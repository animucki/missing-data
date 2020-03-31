#Fit the class model of Lin, McCulloch, and Rosenheck
fit.class <- function(d) {
  key <- as.integer(d$sample[1])
  set.seed(4114L + key)
  flog.debug(paste0('Fitting class model to sample ', key))

  nSubjects <- length(unique(d$subject))
  nTimePoints <- d %>% group_by(subject) %>% summarize(n=n()) %>% pull(n) %>% max
  nClasses <- 3
  dObs <- d %>% filter(r==1)
  
  #Fit ignorable model to find initial values for parameters
  
  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  pars <- list(beta = fixef(m),
               alpha2 = 1e-4, # time parameter of the hazard function - the likelihood has a removable discontinuity at gamma[1]==0 so don't start too close to it
               alpha3 = 0, # treatment parameter
               sigma.b = as.data.frame(VarCorr(m))$sdcor[1],
               sigma = sigma(m),
               theta = 1e-2, # variance of the frailty distribution
               mu = rep(0, nClasses-1),
               eta = rep(0, nClasses-1),
               lambda = rep(1, nClasses) # class-specific baseline hazard
  )

  # smart initialization to make separation of class intercepts more likely
  pars$mu <- c(-1,1) * pars$sigma.b

  previousPars <- Inf
  currentPars <- unlist(pars)
  
  iter <- 1
  
  mcSamples <- 5
  
  minusTwoLogLikelihood <- NA

  dPredictedList <- list()

  nTimesCriterionMet <- 0
  crit <- Inf

  dList <- d %>% group_split(subject)
  doList <- d %>% filter(r==1) %>% group_split(subject)

  while (nTimesCriterionMet < 3 && iter <= 100) {
    
    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ', MC samples: ', mcSamples,', pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ','),
                      ', crit = ', format(crit, digits = 4)) )

    # Monte Carlo Expectation step
    # Initialize the list of simulated datasets
    for (mc in 1:mcSamples) {
      dPredictedList[[mc]] <- matrix(NA_integer_, nrow = nSubjects, ncol = nClasses)
    }

    #simulate
    for (i in 1:nSubjects) {

      # find conditional distribution of ci for this subject
      pr <- rep(NA_real_, nClasses)
      for (k in 1:nClasses) {
        c1 <- pars$lambda[k] * exp(pars$alpha2 * dList[[i]]$time + pars$alpha3 * dList[[i]]$treatment)
        c2 <- - pars$lambda[k] / pars$alpha2 * (exp(pars$alpha2 * dList[[i]]$time)-1) * exp(pars$alpha3 * dList[[i]]$treatment)

        cik <- rep(0, nClasses)
        cik[k] <- 1

        pr[k] <- dmultinormal(x = doList[[i]]$y,
                              mean = pars$beta[1] + pars$beta[2] * doList[[i]]$time + pars$beta[3] * doList[[i]]$treatment + as.vector(c(0,pars$mu) %*% cik),
                              sigma = as.vector( pars$sigma^2 * diag(nrow(doList[[i]])) + pars$sigma.b^2)) *
          prod( (c1/(1-pars$theta * c2))^(1-dList[[i]]$r) * (1-pars$theta * c2)^(-1/pars$theta) ) *
          dmultinomial( x = cik,
                        size = 1,
                        prob = softmax(c(0,pars$eta)))

      }

      ci <- rmultinomial(n = mcSamples, size = 1, prob = pr)

      # Save values
      for (mc in 1:mcSamples) {
        dPredictedList[[mc]][i,] <- ci[mc,]
      }
    }

    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ' E step completed'))
    
    Xtemp <- as.matrix(d %>% select(time, treatment))
    
    # Minimization step
    minusTwoLogLikelihood <- function(x) {
      
      beta <- x[1:3]
      alpha2 <- x[4]
      alpha3 <- x[5]
      sigma.b <- x[6]
      sigma <- x[7]
      theta <- x[8]
      mu <- x[8 + 1:(nClasses-1)]
      eta <- x[8 + (nClasses-1) + 1:(nClasses-1)]
      lambda <- x[8 + 2*(nClasses-1) + 1:nClasses]

      tmp <- lapply(dPredictedList, 
                    function(dObj) {
                      ll <- 0

                      #maybe redo, without the for loop?
                      for (i in 1:nSubjects) {
                        c1 <- lambda[k] * exp(alpha2 * dList[[i]]$time + alpha3 * dList[[i]]$treatment)
                        c2 <- - lambda[k] / alpha2 * (exp(alpha2 * dList[[i]]$time)-1) * exp(alpha3 * dList[[i]]$treatment)

                        ll <- ll +
                          dmultinormal(x = doList[[i]]$y,
                                       mean = beta[1] + beta[2] * doList[[i]]$time + beta[3] * doList[[i]]$treatment + as.vector(c(0,mu) %*% dObj[i,]),
                                       sigma = as.vector( sigma^2 * diag(nrow(doList[[i]])) + sigma.b^2 ),
                                       log = T) +
                          sum( (1-dList[[i]]$r) * log( c1/(1-theta*c2) ) - (1/theta)*log(1 - theta*c2) ) +
                          sum(dmultinomial( x = dObj[i,],
                                            size = 1,
                                            prob = softmax(c(0,eta)),
                                            log = T))
                      }

                      return(ll)
                    })
      
      out <- -2*mean(unlist(tmp))
      
      out
      
    }
    
    minusTwoScore <- function(x) {

      beta <- x[1:3]
      alpha2 <- x[4]
      alpha3 <- x[5]
      sigma.b <- x[6]
      sigma <- x[7]
      theta <- x[8]
      mu <- x[8 + 1:(nClasses-1)]
      eta <- x[8 + (nClasses-1) + 1:(nClasses-1)]
      lambda <- x[8 + 2*(nClasses-1) + 1:nClasses]
      
      out <- lapply(dPredictedList, 
                    function(dObj) {
                      #use mixed grad thing. what was that paper called again?

                      #for missingness and eta, we have very straightforward results. later.

                      grad <- NULL

                      return(grad)
                    }) %>% bind_rows %>% summarize_all(mean) %>% unlist
      -2*out
    }
    
    lowerBounds <- rep(-Inf, length(unlist(pars)))
    lowerBounds[c(6,7,8, 8 + 2*(nClasses-1) + 1:nClasses, length(unlist(pars)))] <- 1e-4 #the three variances, and the lambdas, have to be positive
    
    res <- optim(par=unlist(pars, use.names = F), 
                 fn=minusTwoLogLikelihood,
                 # gr=minusTwoScore,
                 method = 'L-BFGS-B',
                 control = list(trace=3, REPORT=1),
                 lower = lowerBounds, 
                 upper = Inf,
                 hessian = FALSE
    )


    stop()

    if(res$convergence > 0 & iter > 1) flog.error(paste('class likelihood did not converge, code',res$convergence))
    
    pars <- list(beta = res$par[1:3],
                 sigma.b = res$par[4],
                 sigma = res$par[5],
                 mu = res$par[5 + 1:(nClasses-1)],
                 eta = res$par[5 + (nClasses-1) + 1:(nClasses-1)],
                 gamma = res$par[5 + 2*(nClasses-1) + 1:2],
                 llambda = res$par[5 + 2*(nClasses-1) + 2 + 1:nClasses],
                 theta = res$par[5 + 3*(nClasses-1) + 4]
    )
    
    previousMinus2LL <- currentMinus2LL
    previousPars <- currentPars
    currentMinus2LL <- res$value
    currentPars <- unlist(pars)
    
    if(iter < 5)
    {
      mcSamples <- iter
    } else {
      mcSamples <- as.integer(min(mcSamples * 1.5, max(1e7 / nrow(d), 250))) #increase the sample size slowly
    }
    iter <- iter + 1
  }
  
  flog.trace(paste0('Sample ', key, ': EM result for class: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ','), ' , deviance = ', format(currentMinus2LL, digits=7) ) )
  
  x0 <- unlist(pars)
  # hh <- hessian(function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5], x0[9+ 1:(2*nClasses-2)])), x0[c(1,2,3,8,9)])
  hh <- optimHess(par = x0[c(1,2,3,4,5)],
                  fn = function(x) minusTwoLogLikelihood(c(x, x0[5 + 1:(3*(nClasses-1) - 1)])),
                  gr = function(x) minusTwoScore(c(x, x0[5 + 1:(3*(nClasses-1) - 1)]))[c(1,2,3,4,5)]
  )
  
  out <- c(key, pars$beta, pars$sigma.b, pars$sigma, 2*diag(solve(hh)) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}
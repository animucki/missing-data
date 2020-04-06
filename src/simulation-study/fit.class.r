#Fit the class model of Lin, McCulloch, and Rosenheck
fit.class <- function(d) {
  key <- as.integer(d$sample[1])
  set.seed(4114L + key)
  flog.debug(paste0('Fitting class model to sample ', key))

  # Read partial results
  dPartial <- read.csv('/root/balinPartial.csv')
  rowPartial <- 0
  if(key %in% dPartial[,1]) {
    rowPartial <- which(key == dPartial[,1])
  }

  nSubjects <- length(unique(d$subject))
  nTimePoints <- d %>% group_by(subject) %>% summarize(n=n()) %>% pull(n) %>% max
  nClasses <- 3

  # prepare grad0 (initial gradient)
  grad0 <- data.frame(grad.beta1 = 0, grad.beta2 = 0, grad.beta3 = 0,
                      grad.alpha2 = 0, grad.alpha3 = 0,
                      grad.sigma.b = 0, grad.sigma = 0, grad.theta = 0)
  for (k in 2:nClasses) {
    grad0[[paste0('grad.mu',k)]] <- 0
  }
  for (k in 2:nClasses) {
    grad0[[paste0('grad.eta',k)]] <- 0
  }
  for (k in 1:nClasses) {
    grad0[[paste0('grad.lambda',k)]] <- 0
  }
  
  #Fit ignorable model to find initial values for parameters
  
  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  pars <- list(beta = fixef(m),
               alpha2 = 1e-2, # time parameter of the hazard function - the likelihood has a removable discontinuity at alpha2==0 so don't start too close to it
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
    if(rowPartial > 0) {
      flog.trace(paste0('Sample ',key, ': resuming EM from partial result'))
      t <- unlist(dPartial[rowPartial, -1], use.names = F)
      pars <- list(beta = t[1:3],
                 alpha2 = t[4],
                 alpha3 = t[5],
                 sigma.b = t[6],
                 sigma = t[7],
                 theta = t[8],
                 mu = t[8 + 1:(nClasses-1)],
                 eta = t[8 + (nClasses-1) + 1:(nClasses-1)],
                 lambda = t[8 + 2*(nClasses-1) + 1:nClasses])
      mcSamples <- 100
    }

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

      c1wlk <- exp(alpha2 * d$time + alpha3 * d$treatment)
      c2wlk <- - 1/ alpha2 * (exp(alpha2 * d$time)-1) * exp(alpha3 * d$treatment)

      tmp <- lapply(dPredictedList, 
                    function(dObj) {
                      lambdak <- rep(as.vector(dObj %*% lambda), each = nTimePoints)
                      c1 <- c1wlk * lambdak
                      c2 <- c2wlk * lambdak

                      ll <- 0
                      for (i in 1:nSubjects) {
                        ll <- ll +
                          dmultinormal(x = doList[[i]]$y,
                                       mean = beta[1] + beta[2] * doList[[i]]$time + beta[3] * doList[[i]]$treatment + as.vector(c(0,mu) %*% dObj[i,]),
                                       sigma = as.vector( sigma^2 * diag(nrow(doList[[i]])) + sigma.b^2 ),
                                       log = T)
                      }
                          ll <- ll + sum( ((1-d$r) * log( c1/(1-theta*c2) ) - (1/theta)*log(1 - theta*c2))[-seq(from=1, to=nrow(d), by=nTimePoints)] ) +
                          sum(dmultinomial( x = dObj,
                                            size = 1,
                                            prob = softmax(c(0,eta)),
                                            log = T))

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

      # c1wlk <- exp(alpha2 * d$time + alpha3 * d$treatment)
      # c2wlk <- - 1/ alpha2 * (exp(alpha2 * d$time)-1) * exp(alpha3 * d$treatment)

      out <- lapply(dPredictedList, 
                    function(dObj) {

                      # c1 <- c1wlk * lambdak
                      # c2 <- c2wlk * lambdak

                      grad <- grad0
                      for (i in 1:nSubjects) {
                        #precalculate quantities
                        ni <- nrow(doList[[i]])
                        Xi <- cbind(
                          matrix(c(rep(1,ni),
                                   doList[[i]]$time,
                                   doList[[i]]$treatment),
                                 nrow = ni),
                          matrix(rep(dObj[i,-1],
                                     each = ni),
                                 nrow = ni))

                        invSigma <- solve( sigma^2 * diag(ni) + sigma.b^2 )
                        ei <- doList[[i]]$y - Xi %*% c(beta, mu)

                        tau <- invSigma %*% (ei %*% t(ei) - (sigma^2 * diag(ni) + sigma.b^2)) %*% invSigma

                        # grad.beta and grad.mu can be calculated at once (score equations as derived by Jennrich and Schluchter)
                        grad[,c(1,2,3,8 + 1:(nClasses-1))] <- grad[,c(1,2,3,8 + 1:(nClasses-1))] + t(Xi) %*% invSigma %*% ei
                        # then, we have the variance components

                        grad[,6] <- grad[,6] + 1/2 * tr(tau %*% (matrix(2*sigma.b, nrow = ni, ncol = ni)))
                        grad[,7] <- grad[,7] + 1/2 * tr(tau %*% (2*sigma * diag(ni)))

                      }

                      #remaining components can be calculated globally:
                      lambdak <- rep(as.vector(dObj %*% lambda), each = nTimePoints)

                      #grad.alpha2
                      grad[4] <- -sum(
                        ((exp(d$treatment*alpha3) * (exp(d$time*alpha2) - 1) * ((d$r-1)*theta - 1) * lambdak +
                          d$time*alpha2 * ((d$r-1)*alpha2 + exp(d$treatment*alpha3) * (exp(d$time*alpha2) + theta - d$r*theta)*lambdak))/
                          (alpha2*(alpha2 + exp(d$treatment*alpha3)*(exp(d$time*alpha2) - 1)*theta*lambdak))
                        )[-seq(from = 1, to = nrow(d), by = nTimePoints)])

                      #grad.alpha3
                      grad[5] <- -sum(
                        ((d$treatment*((d$r-1)*alpha2 + exp(d$treatment*alpha3) * (exp(d$time*alpha2) - 1)*lambdak))/
                          (alpha2 + exp(d$treatment*alpha3) * (exp(d$time*alpha2) - 1) * theta * lambdak)
                        )[-seq(from = 1, to = nrow(d), by = nTimePoints)])

                      #grad.theta
                      grad[8] <- sum(
                        ((exp(d$treatment*alpha3) * (exp(d$time*alpha2) - 1) * theta * ((d$r-1)*theta -1)*lambdak) /
                          (alpha2 + exp(d$treatment*alpha3) * (exp(d$time*alpha2) - 1) * theta * lambdak) +
                          log(1 + (exp(d$treatment*alpha3) * (exp(d$time*alpha2) - 1) * theta * lambdak)/alpha2)
                        )[-seq(from = 1, to = nrow(d), by = nTimePoints)])/theta^2

                      #grad.eta
                      grad[,8 + (nClasses-1) + 1:(nClasses-1)] <- colSums(dObj[,-1] - rep(softmax(c(0,eta))[-1], each = nrow(dObj)))
                      # print(dObj[,-1] - rep(softmax(c(0,eta))[-1], each = nrow(dObj)))
                      # stop()

                      #grad.lambda
                      dlambdak <- ((1-d$r)*alpha2 - exp(d$treatment *alpha3)*(exp(d$time *alpha2)-1)*lambdak ) /
                        (lambdak * (alpha2 + exp(d$treatment *alpha3)*(exp(d$time *alpha2)-1)*theta*lambdak ))
                      dlambdak <- dlambdak[-seq(from = 1, to = nrow(d), by = nTimePoints)]
                      grad[,8 + 2*(nClasses-1) + 1:nClasses] <-
                        colSums(matrix(rep(dlambdak, times = nClasses), ncol = nClasses) *
                                  dObj[rep(seq_len(nrow(dObj)), each=nTimePoints-1),])



                      return(grad)
                    }) %>% bind_rows %>% summarize_all(mean) %>% unlist

      -2*out
    }
    
    lowerBounds <- rep(-Inf, length(unlist(pars)))
    lowerBounds[c(6,7,8, 8 + 2*(nClasses-1) + 1:nClasses, length(unlist(pars)))] <- 1e-4 #the three variances, and the lambdas, have to be positive

    # ctrl <- list(trace=1, REPORT=10, maxit=max(10*iter,100))
    ctrl <- list(maxit=max(10*iter,100))
    if(nTimesCriterionMet == 2) ctrl$maxit <- 250

    res <- optim(par=unlist(pars, use.names = F), 
                 fn=minusTwoLogLikelihood,
                 gr=minusTwoScore,
                 method = 'L-BFGS-B',
                 control = ctrl,
                 lower = lowerBounds, 
                 upper = Inf,
                 hessian = FALSE
    )

    if(res$convergence == 1 && nTimesCriterionMet >= 2) flog.error('class likelihood did not converge within 250 iterations in final step')
    if(res$convergence > 1) flog.error(paste('class likelihood did not converge, code',res$convergence))
    
    pars <- list(beta = res$par[1:3],
                 alpha2 = res$par[4],
                 alpha3 = res$par[5],
                 sigma.b = res$par[6],
                 sigma = res$par[7],
                 theta = res$par[8],
                 mu = res$par[8 + 1:(nClasses-1)],
                 eta = res$par[8 + (nClasses-1) + 1:(nClasses-1)],
                 lambda = res$par[8 + 2*(nClasses-1) + 1:nClasses]
    )

    previousPars <- currentPars
    currentPars <- unlist(pars)

    mcSamples <- floor(min(mcSamples * 1.2 + 2, 500)) #increase the sample size slowly
    iter <- iter + 1

    #stopping criteria calculation
    crit <- coalesce(mean( (previousPars-currentPars)^2 )/(mean(previousPars^2)+1e-3), Inf)
    if(crit<5e-4) {
      nTimesCriterionMet <- nTimesCriterionMet + 1
    } else {
      nTimesCriterionMet <- 0
    }
    
    #do not repeat if this was a resumed run
    if(rowPartial > 0) break
  }
  
  flog.trace(paste0('Sample ', key, ': EM result for class: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ',')))
  
  x0 <- unlist(pars)
  # hh <- hessian(function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5], x0[9+ 1:(2*nClasses-2)])), x0[c(1,2,3,8,9)])
  hh <- optimHess(par = x0[c(1,2,3,6,7)],
                  fn = function(x) minusTwoLogLikelihood(c(x[c(1,2,3)], x0[c(4,5)], x[c(4,5)], x0[8:length(x0)])),
                  gr = function(x) minusTwoScore(c(x[c(1,2,3)], x0[c(4,5)], x[c(4,5)], x0[8:length(x0)]))[c(1,2,3,6,7)]
  )
  out <- c(key, pars$beta, pars$sigma.b, pars$sigma, sqrt(diag(2*solve(hh))) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')

  flog.trace(paste0('Sample ', key,': full output: ', paste(out, collapse = " ") ))

  as.data.frame(as.list(out))
}

#Fit the class model of Lin, McCulloch, and Rosenheck, but with Bernoulli missingness
fit.hybrid <- function(d) {
  key <- as.integer(d$sample[1])
  set.seed(4110L + key)
  flog.debug(paste0('Fitting hybrid model to sample ', key))

  nSubjects <- length(unique(d$subject))
  nTimePoints <- d %>% group_by(subject) %>% summarize(n=n()) %>% pull(n) %>% max
  nClasses <- 3
  
  #Fit ignorable model to find initial values for parameters
  
  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  pars <- list(beta = fixef(m),
               alpha = c(0,0,0),
               gamma = 0,
               delta = 0,
               sigma.b = as.data.frame(VarCorr(m))$sdcor[1],
               sigma = sigma(m),
               mu = c(-2,3),
               eta = rep(0, nClasses-1)
  )
  
  ##Loop for MCEM
  currentPars <- unlist(pars)
  
  iter <- 1
  
  mcSamples <- 5
  
  minusTwoLogLikelihood <- NA

  dPredictedList <- list()

  nTimesCriterionMet <- 0
  crit <- Inf
  
  while (nTimesCriterionMet < 3 && iter <= 100) {
    
    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ', MC samples: ', mcSamples,', pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ','),
                      ', crit = ', format(crit, digits = 4)) )
    
    # Monte Carlo Expectation step
    # Initialize the list of simulated datasets
    for (mc in 1:mcSamples) {
      dPredictedList[[mc]] <- list(cDrawMat=matrix(NA_integer_,
                                                   nrow = nSubjects,
                                                   ncol = nClasses),
                                   data=d)
      dPredictedList[[mc]]$data$bDraw <- NA_real_
    }

    #simulate
    for (i in 1:nSubjects) {
      di <- d %>% filter(subject==i)

      #simulate bi first
      bVec <- arms(n_samples = mcSamples, lower = -10, upper = 10, metropolis = TRUE,
                   log_pdf = function(bi) {
                     out <- 0
                     for (k in 1:nClasses){
                       cik <- rep.int(0,nClasses)
                       cik[k] <- 1

                       out <- out +
                         prod(dnorm( #outcome
                           x = di$y,
                           mean = pars$beta[1] + di$time * pars$beta[2] + di$treatment * pars$beta[3] +
                             bi + c(0, pars$mu)[k],
                           sd = pars$sigma), na.rm = T) *
                           prod(dbinom( #missingness
                             x = di$r,
                             size = 1,
                             prob = plogis(pars$alpha[1] + di$time * pars$alpha[2] + di$treatment * pars$alpha[3] +
                                             pars$gamma * bi + pars$delta * c(0, pars$mu)[k])
                           )[-1]) *
                           dnorm( #(normal) random intercept
                             x = bi,
                             sd = pars$sigma.b) *
                           dmultinomial( #class-specific intercept
                             x = cik,
                             size = 1,
                             prob = softmax(c(0, pars$eta)))
                     }

                     return(log(out))
                   })

      for (mc in 1:mcSamples) {
        # Find the conditional class membership distribution
        ciConditional <- rep(NA_real_, nClasses)
        for (k in 1:nClasses) {
          cik <- rep.int(0,nClasses)
          cik[k] <- 1

          ciConditional[k] <-
            prod(dnorm( #outcome
              x = di$y,
              mean = pars$beta[1] + di$time * pars$beta[2] + di$treatment * pars$beta[3] +
                bVec[mc] + c(0, pars$mu)[k],
              sd = pars$sigma), na.rm = T) *
              prod(dbinom( #missingness
                x = di$r,
                size = 1,
                prob = plogis(pars$alpha[1] + di$time * pars$alpha[2] + di$treatment * pars$alpha[3] +
                                pars$gamma * bVec[mc] + pars$delta * c(0, pars$mu)[k])
              )[-1]) *
              dnorm( #(normal) random intercept
                x = bVec[mc],
                sd = pars$sigma.b) *
              dmultinomial( #class-specific intercept
                x = cik,
                size = 1,
                prob = softmax(c(0, pars$eta)))

        }

        # Simulate class and save values
        dPredictedList[[mc]]$cDrawMat[i,] <- rmultinomial(n=1, size=1, prob=ciConditional)
        dPredictedList[[mc]]$data$bDraw[d$subject==i] <- bVec[mc]
      }
    }

    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ' E step completed'))
    
    Xtemp <- as.matrix(d %>% select(time, treatment))
    
    # Minimization stepPredictedList
    minusTwoLogLikelihood <- function(x) {
      
      if(length(x) != 10 + 2*(nClasses-1)) flog.error('Incorrect input in logLikelihood for hybrid model')

      beta <- x[1:3]
      alpha <- x[4:6]
      gamma <- x[7]
      delta <- x[8]
      sigma.b <- x[9]
      sigma <- x[10]
      mu <- x[10 + 1:(nClasses-1)]
      eta <- x[10 + (nClasses-1) + 1:(nClasses-1)]

      tmp <- lapply(dPredictedList, 
                    function(dObj) {
                      out <- sum(dnorm( #outcome
                        x = dObj$data$y,
                        mean = beta[1] + Xtemp %*% beta[c(2,3)] + dObj$data$bDraw + rep(dObj$cDrawMat %*% c(0, mu), each = nTimePoints),
                        sd = sigma,
                        log = T), na.rm = TRUE) +
                        sum(dbinom( #missingness
                          x= dObj$data$r,
                          size = 1,
                          prob = plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] + gamma * dObj$data$bDraw + delta * rep(dObj$cDrawMat %*% c(0, mu), each = nTimePoints)),
                          log = T)[-seq(from=1, to=nrow(dObj$data), by=nTimePoints)])+
                        #second, add up subject-level contributions
                        sum(dnorm( #(normal) random intercept
                          x = dObj$data$bDraw[seq(1, nrow(dObj$data), by = nTimePoints)],
                          sd = sigma.b,
                          log = T)) +
                        sum(dmultinomial( #class-specific intercept
                          x = dObj$cDrawMat,
                          size = 1,
                          prob = softmax(c(0,eta)),
                          log = T))
                      out
                    })
      
      out <- -2*mean(unlist(tmp))

      out
      
    }
    
    minusTwoScore <- function(x) {
      
      if(length(x) != 10 + 2*(nClasses-1)) flog.error('Incorrect input in logLikelihood for parametric model')

      beta <- x[1:3]
      alpha <- x[4:6]
      gamma <- x[7]
      delta <- x[8]
      sigma.b <- x[9]
      sigma <- x[10]
      mu <- x[10 + 1:(nClasses-1)]
      eta <- x[10 + (nClasses-1) + 1:(nClasses-1)]
      
      out <- lapply(dPredictedList, 
                    function(dObj) {
                      mcVec <- rep(dObj$cDrawMat %*% c(0, mu), each = nTimePoints)
                      normalResidual <- dObj$data$y - beta[1] - Xtemp %*% beta[c(2,3)] - dObj$data$bDraw - mcVec
                      bernoulliResidual <- dObj$data$r - plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] +
                                                                  gamma * dObj$data$bDraw + delta * mcVec)

                      grad <- data.frame(
                         grad.beta1 = sum( normalResidual, na.rm = T )/sigma^2,
                         grad.beta2 = sum( dObj$data$time * normalResidual, na.rm = T )/sigma^2,
                         grad.beta3 = sum( dObj$data$treatment * normalResidual, na.rm = T )/sigma^2,
                         grad.alpha1 = sum( (bernoulliResidual)[-seq(from=1, to=nrow(dObj$data), by=nTimePoints)] ),
                         grad.alpha2 = sum( (dObj$data$time * bernoulliResidual)[-seq(from=1, to=nrow(dObj$data), by=nTimePoints)] ),
                         grad.alpha3 = sum( (dObj$data$treatment * bernoulliResidual)[-seq(from=1, to=nrow(dObj$data), by=nTimePoints)] ),
                         grad.gamma = sum( (dObj$data$bDraw * bernoulliResidual)[-seq(from=1, to=nrow(dObj$data), by=nTimePoints)] ),
                         grad.delta = sum( (mcVec * bernoulliResidual)[-seq(from=1, to=nrow(dObj$data), by=nTimePoints)] ),
                         grad.sigma.b = sum( (dObj$data$bDraw[seq(1, nrow(dObj$data), by = nTimePoints)]^2 - sigma.b^2 )/sigma.b^3),
                         grad.sigma = sum( normalResidual^2 - sigma^2, na.rm = T )/sigma^3)
                      
                      #two separate loops ensure the correct order in the output!
                      for(k in 2:nClasses) {
                        grad[[paste0('grad.mu',k)]] <- sum( normalResidual * rep(dObj$cDrawMat[,k], each = nTimePoints), na.rm = T )/sigma^2 +
                          sum( (delta * rep(dObj$cDrawMat[,k], each = nTimePoints) * bernoulliResidual)[-seq(from=1, to=nrow(dObj$data), by=nTimePoints)] )
                      }
                      
                      for(k in 2:nClasses) {
                        grad[[paste0('grad.eta',k)]] <- sum( dObj$cDrawMat[,k] - softmax(c(0,eta))[k] )
                      }

                      return(grad)
                    }) %>% bind_rows %>% summarize_all(mean) %>% unlist

      # print(c(value=-2*out['grad.delta'], relerr=(-2*out/grad(minusTwoLogLikelihood,x))['grad.delta']))

      -2*out

      
    }

    res <- optim(par=unlist(pars, use.names = F), 
                 fn=minusTwoLogLikelihood,
                 gr=minusTwoScore,
                 method = 'L-BFGS-B',
                 # control = list(trace=1, REPORT=1),
                 lower = c(rep(-Inf, 8), 1e-4, 1e-4, rep(-Inf, 2*(nClasses-1))),
                 upper = Inf,
                 hessian = FALSE
    )
    
    if(res$convergence > 0) flog.warn(paste('SPM+class likelihood did not converge, code',res$convergence))

    pars <- list(beta = res$par[1:3],
                 alpha = res$par[4:6],
                 gamma = res$par[7],
                 delta = res$par[8],
                 sigma.b = res$par[9],
                 sigma = res$par[10],
                 mu = res$par[10 + 1:(nClasses-1)],
                 eta = res$par[10 + (nClasses-1) + 1:(nClasses-1)])

    previousPars <- currentPars
    currentPars <- unlist(pars)

    mcSamples <- as.integer(min(mcSamples * 1.25, 250)) #increase the sample size slowly
    iter <- iter + 1

    #stopping criteria calculation
    crit <- coalesce(mean( (previousPars-currentPars)^2 )/(mean(previousPars^2)+1e-3), Inf)
    if(crit<1e-4) {
      nTimesCriterionMet <- nTimesCriterionMet + 1
    } else {
      nTimesCriterionMet <- 0
    }
  }
  
  flog.trace(paste0('Sample ', key, ': EM result for spm+class: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ',') ) )
  
  x0 <- unlist(pars)
  hh <- optimHess(par = x0[c(1,2,3,9,10)],
                  fn = function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:8], x[4:5], x0[10+ 1:(2*nClasses-2)])),
                  gr = function(x) minusTwoScore(c(x[1:3], x0[4:8], x[4:5], x0[10+ 1:(2*nClasses-2)]))[c(1,2,3,9,10)]
                  )

  out <- c(key, pars$beta, pars$sigma.b, pars$sigma, sqrt(diag(2*solve(hh))) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}
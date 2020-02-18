fit.class <- function(d) {
  key <- as.integer(d$sample[1])
  set.seed(4110L + key)
  flog.debug(paste0('Fitting class model to sample ', key))
  
  nTimePoints <- d %>% group_by(subject) %>% summarize(n=n()) %>% pull(n) %>% max
  nClasses <- 3
  
  #Fit ignorable model to find initial values for parameters
  
  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  pars <- list(beta = fixef(m),
               sigma.b = as.data.frame(VarCorr(m))$sdcor[1],
               sigma = sigma(m),
               mu = rep(0, nClasses-1),
               eta = rep(0, nClasses-1),
               gamma = c(0,1e-4), #parameters of the hazard function - the likelihood has a removable discontinuity at gamma[2]==0 so don't start too close to it
               llambda = rep(0, nClasses), #log(baseline hazard)
               theta = 1 #variance of the frailty distribution
  )
  
  #calculate tSince (time-to-event where the event is missingness)
  d <- d %>% mutate(tEvent = case_when(
    r == 0 ~ time,
    r == 1 ~ NA_real_
  ))
  
  ##Loop for stochastic-EM
  previousMinus2LL <- Inf
  currentMinus2LL <- Inf
  
  previousPars <- Inf
  currentPars <- unlist(pars)
  
  iter <- 1
  
  mcSamples <- 8
  dPredictedList <- list()
  Xtemp <- NA
  
  minusTwoLogLikelihood <- NA
  
  while (coalesce(abs(previousMinus2LL-currentMinus2LL), Inf) > 1e-6 && 
         coalesce(mean( (previousPars-currentPars)^2 ), Inf) > 1e-4 && 
         iter <= 100) {
    
    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ', MC samples: ', mcSamples,', pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ','),', normChange = ', format(sum( (previousPars-currentPars)^2 ), digits = 4)) )
    
    # Stochastic step
    for (mc in 1:mcSamples) {
      #cDrawMat, then assign to list object
      dPredictedList[[mc]] <- list(cDrawMat=NULL, data=NULL)
      dPredictedList[[mc]]$cDrawMat <- rmultinomial(n=length(unique(d$subject)), size=1, prob = softmax(c(0, pars$eta)))
      dPredictedList[[mc]]$data <- d %>% mutate(
        bDraw = rep(rnorm(n=length(unique(d$subject)), sd=pars$sigma.b), each = nTimePoints),
        mcDraw = rep(dPredictedList[[mc]]$cDrawMat %*% c(0, pars$mu), each = nTimePoints),
        eDraw = rnorm(n=nrow(d), sd=pars$sigma),
        yPred = case_when(
          r == 1 ~ y,
          r == 0 ~ bDraw + mcDraw + pars$beta[1] + pars$beta[2] * time + pars$beta[3] * treatment + eDraw
        ),
        omegaDraw = rep(rgamma(n=length(unique(d$subject)), shape=1/pars$theta, scale=pars$theta), each=nTimePoints))
    }
    
    Xtemp <- as.matrix(d %>% select(time, treatment))
    
    # Minimization step
    minusTwoLogLikelihood <- function(x) {
      
      beta <- x[1:3]
      sigma.b <- x[4]
      sigma <- x[5]
      mu <- x[5 + 1:(nClasses-1)]
      eta <- x[5 + (nClasses-1) + 1:(nClasses-1)]
      gamma <- x[5 + 2*(nClasses-1) + 1:2]
      llambda <- x[5 + 2*(nClasses-1) + 2 + 1:nClasses]
      theta <- x[5 + 3*(nClasses-1) + 4]
      
      tmp <- lapply(dPredictedList, 
                    function(dObj) {
                      #first, add up observation-level loglikelihood contributions
                      ll <- sum(
                        dnorm( #outcome
                          x = dObj$data$yPred, 
                          mean = beta[1] + Xtemp %*% beta[c(2,3)] + dObj$data$bDraw + rep(dObj$cDrawMat %*% c(0, mu), each = nTimePoints), 
                          sd = sigma, 
                          log = T) 
                      ) +
                        #second, add up subject-level contributions 
                        sum(
                          dnorm( #(normal) random intercept
                            x = dObj$data$bDraw[seq(1, nrow(dObj$data), by = nTimePoints)], 
                            sd = sigma.b, 
                            log = T) +
                            dmultinomial( #class-specific intercept 
                              x = dObj$cDrawMat,
                              size = 1,
                              prob = softmax(c(0,eta)),
                              log = T) +
                            dgamma( #frailty
                              x = dObj$data$omegaDraw[seq(1, nrow(dObj$data), by = nTimePoints)],
                              shape = 1/theta,
                              scale = theta,
                              log = T
                            )
                        )
                      
                      #add contribution from missingness
                      ll <- ll + sum(
                        rep(dObj$cDrawMat %*% llambda, each = nTimePoints) + log(dObj$data$omegaDraw) + dObj$data$treatment * gamma[1] + dObj$data$tEvent * gamma[2] -
                          exp(dObj$data$treatment * gamma[1]) * (-1+exp(dObj$data$tEvent * gamma[2]))/gamma[2] *
                          rep(dObj$cDrawMat %*% exp(llambda), each = nTimePoints) * dObj$data$omegaDraw,
                        na.rm = T)
                      
                      return(ll)
                    })
      
      out <- -2*mean(unlist(tmp))
      
      out
      
    }
    
    minusTwoScore <- function(x) {
      
      beta <- x[1:3]
      sigma.b <- x[4]
      sigma <- x[5]
      mu <- x[5 + 1:(nClasses-1)]
      eta <- x[5 + (nClasses-1) + 1:(nClasses-1)]
      gamma <- x[5 + 2*(nClasses-1) + 1:2]
      llambda <- x[5 + 2*(nClasses-1) + 2 + 1:nClasses]
      theta <- x[5 + 3*(nClasses-1) + 4]
      
      out <- lapply(dPredictedList, 
                    function(dObj) {
                      #first, add up observation-level loglikelihood contributions
                      normalResidual <- dObj$data$yPred - beta[1] - Xtemp %*% beta[c(2,3)] - 
                        dObj$data$bDraw - rep(dObj$cDrawMat %*% c(0, mu), each = nTimePoints)
                      
                      grad <- data.frame(
                        grad.beta1 = sum( normalResidual )/sigma^2,
                        grad.beta2 = sum( dObj$data$time * normalResidual )/sigma^2,
                        grad.beta3 = sum( dObj$data$treatment * normalResidual )/sigma^2,
                        grad.sigma.b = sum( (dObj$data$bDraw[seq(1, nrow(dObj$data), by = nTimePoints)]^2 - sigma.b^2 )/sigma.b^3),
                        grad.sigma = sum( normalResidual^2 - sigma^2 )/sigma^3)
                      
                      #the separate loops ensure the correct order in the output!
                      for(k in 2:nClasses) {
                        grad[[paste0('grad.mu',k)]] <- sum( normalResidual * rep(dObj$cDrawMat[,k], each = nTimePoints) )/sigma^2
                      }
                      
                      for(k in 2:nClasses) {
                        grad[[paste0('grad.eta',k)]] <- sum( dObj$cDrawMat[,k] - softmax(c(0,eta))[k] )
                      }
                      
                      grad[['grad.gamma1']] <- sum(
                        dObj$data$treatment * (1 - exp(dObj$data$treatment * gamma[1]) * (-1+exp(dObj$data$tEvent * gamma[2]))/gamma[2] *
                                                 rep(dObj$cDrawMat %*% exp(llambda), each = nTimePoints) * dObj$data$omegaDraw),
                        na.rm = T)
                      grad[['grad.gamma2']] <- sum(
                        dObj$data$time - exp(dObj$data$treatment * gamma[1] ) * (1+ exp(dObj$data$tEvent * gamma[2]) * (-1 + dObj$data$tEvent * gamma[2])) *
                          rep(dObj$cDrawMat %*% exp(llambda), each = nTimePoints) * dObj$data$omegaDraw / gamma[2]^2,
                        na.rm = T)
                      
                      for(k in 1:nClasses) {
                        grad[[paste0('grad.llambda',k)]] <- sum(
                          rep(dObj$cDrawMat[,k], each = nTimePoints)*
                          (1 - rep(dObj$cDrawMat[,k] * exp(llambda[k]), each = nTimePoints) * exp(dObj$data$treatment * gamma[1]) *
                            (-1+exp(dObj$data$tEvent * gamma[2]))/gamma[2] * dObj$data$omegaDraw),
                          na.rm = T)
                      }
                      
                      grad[['grad.theta']] <- sum((-1 + dObj$data$omegaDraw[seq(1, nrow(dObj$data), by = nTimePoints)] +log(theta) - 
                                                 log(dObj$data$omegaDraw[seq(1, nrow(dObj$data), by = nTimePoints)]) + digamma(1/theta))/theta^2)
                      
                      return(grad)
                    }) %>% bind_rows %>% summarize_all(mean) %>% unlist
      -2*out
    }
    
    lowerBounds <- rep(-Inf, length(unlist(pars)))
    lowerBounds[c(4,5,length(unlist(pars)))] <- 1e-4 #the two variances, and theta
    
    res <- optim(par=unlist(pars, use.names = F), 
                 fn=minusTwoLogLikelihood,
                 gr=minusTwoScore,
                 method = 'L-BFGS-B',
                 # control = list(trace=3, REPORT=1),
                 lower = lowerBounds, 
                 upper = Inf,
                 hessian = FALSE
    )
    
    if(res$convergence > 0) flog.error(paste('class likelihood did not converge, code',res$convergence))
    
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
    mcSamples <- as.integer(min(mcSamples * log(10)/log(2), max(1e7 / nrow(d), 250)))
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
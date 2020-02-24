#Fit the class model of Lin, McCulloch, and Rosenheck, but with Bernoulli missingness
fit.hybrid <- function(d) {
  key <- as.integer(d$sample[1])
  set.seed(4110L + key)
  flog.debug(paste0('Fitting hybrid model to sample ', key))
  
  nTimePoints <- d %>% group_by(subject) %>% summarize(n=n()) %>% pull(n) %>% max
  nClasses <- 3
  
  #Fit ignorable model to find initial values for parameters
  
  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  pars <- list(beta = fixef(m),
               sigma.b = as.data.frame(VarCorr(m))$sdcor[1],
               sigma = sigma(m),
               alpha = c(0,0,0),
               gamma = 0,
               delta = 0,
               mu = rep(0, nClasses-1),
               eta = rep(0, nClasses-1)
  )
  
  ##Loop for stochastic-EM
  previousMinus2LL <- Inf
  currentMinus2LL <- Inf
  
  previousPars <- Inf
  currentPars <- unlist(pars)
  
  iter <- 1
  
  mcSamples <- 1
  dPredictedList <- list()
  Xtemp <- NA
  
  minusTwoLogLikelihood <- NA
  
  while (coalesce(abs(previousMinus2LL-currentMinus2LL), Inf) > 1e-6 && 
         coalesce(mean( (previousPars-currentPars)^2 ), Inf) > 1e-5 && 
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
        ))
    }
    
    Xtemp <- as.matrix(d %>% select(time, treatment))
    
    # Minimization step
    minusTwoLogLikelihood <- function(x) {
      
      if(length(x) != 9 + 2*(nClasses-1)) flog.error('Incorrect input in logLikelihood for hybrid model')
      
      beta <- x[1:3]
      sigma.b <- x[4]
      sigma <- x[5]
      alpha <- x[6:8]
      gamma <- x[9]
      delta <- x[10]
      mu <- x[10 + 1:(nClasses-1)]
      eta <- x[10 + (nClasses-1) + 1:(nClasses-1)]
      
      tmp <- lapply(dPredictedList, 
                    function(dObj) {
                      #first, add up observation-level loglikelihood contributions
                      sum(
                        dnorm( #outcome
                          x = dObj$data$yPred, 
                          mean = beta[1] + Xtemp %*% beta[c(2,3)] + dObj$data$bDraw + rep(dObj$cDrawMat %*% c(0, mu), each = nTimePoints), 
                          sd = sigma, 
                          log = T) +
                          dbinom( #missingness
                            x= dObj$data$r, 
                            size = 1, 
                            prob = plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] + gamma * dObj$data$bDraw + delta * rep(dObj$cDrawMat %*% c(0, mu), each = nTimePoints)), 
                            log = T)
                      ) +  
                        #second, add up subject-level contributions ()
                        sum(
                          dnorm( #(normal) random intercept
                            x = dObj$data$bDraw[seq(1, nrow(dObj$data), by = nTimePoints)], 
                            sd = sigma.b, 
                            log = T) +
                            +
                            dmultinomial( #class-specific intercept
                              x = dObj$cDrawMat,
                              size = 1,
                              prob = softmax(c(0,eta)),
                              log = T)
                        )
                    })
      
      out <- -2*mean(unlist(tmp))
      
      out
      
    }
    
    minusTwoScore <- function(x) {
      
      if(length(x) != 9 + 2*(nClasses-1)) flog.error('Incorrect input in logLikelihood for parametric model')
      
      beta <- x[1:3]
      alpha <- x[4:6]
      gamma <- x[7]
      sigma.b <- x[8]
      sigma <- x[9]
      mu <- x[9 + 1:(nClasses-1)]
      eta <- x[9 + (nClasses-1) + 1:(nClasses-1)]
      
      out <- lapply(dPredictedList, 
                    function(dObj) {
                      #first, add up observation-level loglikelihood contributions
                      normalResidual <- dObj$data$yPred - beta[1] - Xtemp %*% beta[c(2,3)] - 
                        dObj$data$bDraw - rep(dObj$cDrawMat %*% c(0, mu), each = nTimePoints)
                      bernoulliResidual <- dObj$data$r - plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] + gamma * dObj$data$bDraw )
                      
                      grad <- data.frame(
                         grad.beta1 = sum( normalResidual )/sigma^2,
                         grad.beta2 = sum( dObj$data$time * normalResidual )/sigma^2,
                         grad.beta3 = sum( dObj$data$treatment * normalResidual )/sigma^2,
                         grad.alpha1 = sum( bernoulliResidual ),
                         grad.alpha2 = sum( dObj$data$time * bernoulliResidual ),
                         grad.alpha3 = sum( dObj$data$treatment * bernoulliResidual ),
                         grad.gamma = sum( dObj$data$bDraw * bernoulliResidual ),
                         grad.sigma.b = sum( (dObj$data$bDraw[seq(1, nrow(dObj$data), by = nTimePoints)]^2 - sigma.b^2 )/sigma.b^3),
                         grad.sigma = sum( normalResidual^2 - sigma^2 )/sigma^3)
                      
                      #two separate loops ensure the correct order in the output!
                      for(k in 2:nClasses) {
                        grad[[paste0('grad.mu',k)]] <- sum( normalResidual * rep(dObj$cDrawMat[,k], each = nTimePoints) )/sigma^2
                      }
                      
                      for(k in 2:nClasses) {
                        grad[[paste0('grad.eta',k)]] <- sum( dObj$cDrawMat[,k] - softmax(c(0,eta))[k] )
                      }
                      return(grad)
                    }) %>% bind_rows %>% summarize_all(mean) %>% unlist
      
      -2*out
      
    }

    res <- optim(par=unlist(pars, use.names = F), 
                 fn=minusTwoLogLikelihood,
                 gr=minusTwoScore,
                 method = 'L-BFGS-B',
                 #control = list(trace=3, REPORT=1),
                 lower = c(rep(-Inf, 7), 1e-4, 1e-4, rep(-Inf, 2*(nClasses-1))),
                 upper = Inf,
                 hessian = FALSE
    )
    
    if(res$convergence > 0) flog.error(paste('SPM+class likelihood did not converge, code',res$convergence))
    
    pars <- list(beta = res$par[1:3],
                 alpha = res$par[4:6],
                 gamma = res$par[7],
                 sigma.b = res$par[8],
                 sigma = res$par[9],
                 mu = res$par[9 + 1:(nClasses-1)],
                 eta = res$par[9 + (nClasses-1) + 1:(nClasses-1)])
    
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
  
  flog.trace(paste0('Sample ', key, ': EM result for spm+class: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ','), ' , deviance = ', format(currentMinus2LL, digits=7) ) )
  
  x0 <- unlist(pars)
  # hh <- hessian(function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5], x0[9+ 1:(2*nClasses-2)])), x0[c(1,2,3,8,9)])
  hh <- optimHess(par = x0[c(1,2,3,8,9)],
                  fn = function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5], x0[9+ 1:(2*nClasses-2)])),
                  gr = function(x) minusTwoScore(c(x[1:3], x0[4:7], x[4:5], x0[9+ 1:(2*nClasses-2)]))[c(1,2,3,8,9)]
                  )

  out <- c(key, pars$beta, pars$sigma.b, pars$sigma, 2*diag(solve(hh)) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}
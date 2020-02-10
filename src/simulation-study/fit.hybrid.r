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
               alpha = c(0,0,0),
               gamma = 0,
               sigma.b = as.data.frame(VarCorr(m))$sdcor[1],
               sigma = sigma(m),
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
         coalesce(sum( (previousPars-currentPars)^2 ), Inf) > 1e-3 && 
         iter <= 100) {
    
    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ', MC samples: ', mcSamples,', pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ','),', normChange = ', format(sum( (previousPars-currentPars)^2 ), digits = 4)) )
    
    # Stochastic step
    #once done: scan for dPredicted usages
    for (mc in 1:mcSamples) {
      dPredictedList[[mc]]  <- d %>% mutate(
        bDraw = rep(rnorm(n=length(unique(d$subject)), sd=pars$sigma.b), each = nTimePoints),
        cDraw = rep(as.integer(rmultinomial(n=length(unique(d$subject)), size=1, prob = softmax(c(0, pars$eta))) %*% 1:nClasses), each = nTimePoints),
        mcDraw = outer(cDraw, 1:nClasses, function(x,y) as.integer(x==y)) %*% c(0, pars$mu),
        eDraw = rnorm(n=nrow(d), sd=pars$sigma),
        yPred = case_when(
          r == 1 ~ y,
          r == 0 ~ bDraw + mcDraw + pars$beta[1] + pars$beta[2] * time + pars$beta[3] * treatment + eDraw
        ))
    }
    
    Xtemp <- as.matrix(d %>% select(time, treatment))
    
    # Minimization step
    minusTwoLogLikelihood <- function(x) {
      
      if(length(x) != 9 + 2*(nClasses-1)) flog.error('Incorrect input in logLikelihood for parametric model')
      
      beta <- x[1:3]
      alpha <- x[4:6]
      gamma <- x[7]
      sigma.b <- x[8]
      sigma <- x[9]
      mu <- x[9 + 1:(nClasses-1)]
      eta <- x[9 + (nClasses-1) + 1:(nClasses-1)]
      
      tmp <- lapply(dPredictedList, 
                    function(dPredicted) {
                      #first, add up observation-level loglikelihood contributions
                      sum(
                        dnorm( #outcome
                          x = dPredicted$yPred, 
                          mean = beta[1] + Xtemp %*% beta[c(2,3)] + dPredicted$bDraw + 
                            outer(dPredicted$cDraw, 1:nClasses, function(x,y) as.integer(x==y)) %*% c(0, mu), 
                          sd = sigma, 
                          log = T) +
                          dbinom( #missingness
                            x= dPredicted$r, 
                            size = 1, 
                            prob = plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] + gamma * dPredicted$bDraw ), 
                            log = T)
                      ) +  
                        #second, add up subject-level contributions ()
                        sum(
                          dnorm( #(normal) random intercept
                            x = dPredicted$bDraw[seq(1, nrow(dPredicted), by = nTimePoints)], 
                            sd = sigma.b, 
                            log = T) +
                            +
                            dmultinomial( #class-specific intercept
                              x = outer(dPredicted$cDraw, 1:nClasses, function(x,y) as.integer(x==y)),
                              size = 1,
                              prob = softmax(c(0,eta)),
                              log = T)
                        )
                    })
      
      out <- -2*mean(unlist(tmp))
      
      # if( !is.finite(out) ) print(unlist(tmp)[!is.finite(unlist(tmp))])
      
      out
      
    }
    
    res <- optim(unlist(pars, use.names = F), minusTwoLogLikelihood,
                 method = 'L-BFGS-B',
                 control = list(trace=1, REPORT=1),
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
    mcSamples <- min(mcSamples * 5, 100)
    iter <- iter + 1
  }
  
  flog.trace(paste0('Sample ', key, ': EM result for spm+class: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ','), ' , deviance = ', format(currentMinus2LL, digits=7) ) )
  
  x0 <- unlist(pars)
  hh <- hessian(function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5], x0[9+ 1:(2*nClasses-2)])), x0[c(1,2,3,8,9)])
  
  out <- c(key, pars$beta, pars$sigma.b, pars$sigma, 2*diag(solve(hh)) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}
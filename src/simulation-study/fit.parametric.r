fit.parametric <- function(d) {
  key <- as.integer(d$sample[1])
  set.seed(1410L + key)
  flog.debug(paste0('Fitting parametric model to sample ', key))
  
  nTimePoints <- d %>% group_by(subject) %>% summarize(n=n()) %>% pull(n) %>% max
  
  #Fit ignorable model to find initial values for parameters
  
  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  pars <- list(beta = fixef(m),
               alpha = c(0,0,0),
               gamma = 0,
               sigma.b = as.data.frame(VarCorr(m))$sdcor[1],
               sigma = sigma(m)
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
        eDraw = rnorm(n=nrow(d), sd=pars$sigma),
        yPred = case_when(
          r == 1 ~ y,
          r == 0 ~ bDraw + pars$beta[1] + pars$beta[2] * time + pars$beta[3] * treatment + eDraw
        ))
    }
    
    Xtemp <- as.matrix(d %>% select(time, treatment))
    
    # Minimization step
    minusTwoLogLikelihood <- function(x) {
      
      if(length(x) != 9) flog.error('Incorrect input in logLikelihood for parametric model')
      
      beta <- x[1:3]
      alpha <- x[4:6]
      gamma <- x[7]
      sigma.b <- x[8]
      sigma <- x[9]
      
      temp <- lapply(dPredictedList, 
                       function(dPredicted) {
                         sum(dnorm(dPredicted$yPred - beta[1] - Xtemp %*% beta[c(2,3)] - dPredicted$bDraw, sd = sigma, log = T)+
                           dbinom( dPredicted$r, size = 1, prob = plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] + gamma * dPredicted$bDraw ), log = T)) +
                           sum(dnorm(dPredicted$bDraw[seq(1, nrow(dPredicted), by = nTimePoints)], sd = sigma.b, log = T))
                       })
      
      out <- -2*mean(unlist(temp))
      
      # if( !is.finite(out) ) print(unlist(temp)[!is.finite(unlist(temp))])
      
      out
      
    }
    
    res <- optim(unlist(pars, use.names = F), minusTwoLogLikelihood,
                 method = 'L-BFGS-B',
                 lower = c(rep(-Inf, 7), 1e-4, 1e-4),
                 upper = Inf,
                 hessian = FALSE
                 )
    
    if(res$convergence > 0) flog.error(paste('Parametric model likelihood did not converge, code',res$convergence))
    
    pars <- list(beta = res$par[1:3],
                 alpha = res$par[4:6],
                 gamma = res$par[7],
                 sigma.b = res$par[8],
                 sigma = res$par[9])
    
    previousMinus2LL <- currentMinus2LL
    previousPars <- currentPars
    currentMinus2LL <- res$value
    currentPars <- unlist(pars)
    mcSamples <- min(mcSamples * 5, 100)
    iter <- iter + 1
  }
  
  flog.trace(paste0('Sample ', key, ': EM result: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ','), ' , deviance = ', format(currentMinus2LL, digits=7) ) )
  
  flog.trace('Computing Hessian...')
  x0 <- unlist(pars)
  hh <- hessian(function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5])), x0[c(1,2,3,8,9)])
  
  flog.debug(paste0('Sample ', key, ' done.'))
  
  out <- c(key, pars$beta, pars$sigma.b, pars$sigma, 2*diag(solve(hh)) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}
#Fit the shared-parameter model
fit.parametric <- function(d) {
  #IMPORTANT: this function assumes that the observations are ordered by time!!!
  
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
  
  while (coalesce(abs(currentMinus2LL-previousMinus2LL)/previousMinus2LL, Inf) > 1e-4 && 
         # coalesce(mean( (previousPars-currentPars)^2 ), Inf) > 1e-5 &&
         iter <= 100) {
    
    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ', MC samples: ', mcSamples,', pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ','),
                      ', normChange = ', format(sum( (previousPars-currentPars)^2 ), digits = 4),
                      ', deviance = ', format(currentMinus2LL, digits=7) ))
    
    # Monte Carlo step
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
                         sum(dmvnorm(
                           x = matrix(dPredicted$yPred - beta[1] - Xtemp %*% beta[c(2,3)], ncol = nTimePoints, byrow = T),
                           sigma = sigma.b^2 + sigma^2 * diag(nTimePoints),
                           log = T
                         ))+
                           sum(dbinom( 
                             x=dPredicted$r, 
                             size = 1, 
                             prob = plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] + gamma * dPredicted$bDraw ),
                             log = T))
                       })
      
      out <- -2*mean(unlist(temp))
      
      out
      
    }
    
    minusTwoScore <- function(x) {
      
      if(length(x) != 9) flog.error('Incorrect input in score for parametric model')
      
      beta <- x[1:3]
      alpha <- x[4:6]
      gamma <- x[7]
      sigma.b <- x[8]
      sigma <- x[9]
      
      temp <- lapply(dPredictedList, 
                     function(dPredicted) {
                       
                       #Pre-calculation of various parts of the score
                       bernoulliResidual <- dPredicted$r - plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] + gamma * dPredicted$bDraw )
                       invSigma <- solve(sigma.b^2 + sigma^2*diag(nTimePoints))
                       e <- matrix(dPredicted$yPred - beta[1] - Xtemp %*% beta[c(2,3)], nrow = nTimePoints)
                       invSigma.e <- as.vector(invSigma %*% e)
                       
                       tau <- matrix(0, nrow = nTimePoints, ncol = nTimePoints)
                       grad.sigma.b <- 0
                       grad.sigma <- 0
                       for(i in 1:ncol(e)) {
                         #Pre-calculation of common component of variance parameter score
                         tau <- invSigma %*% (e[,i] %o% e[,i] - (sigma.b^2 + sigma^2*diag(nTimePoints))) %*% invSigma
                         
                         grad.sigma.b <- grad.sigma.b + sum(diag( tau %*% matrix(sigma.b, nrow=nTimePoints, ncol = nTimePoints)))
                         grad.sigma <- grad.sigma + sum(diag( tau )) * sigma
                       }
                       
                       data.frame(
                         grad.beta1 = sum( invSigma.e ),
                         grad.beta2 = sum( dPredicted$time * invSigma.e),
                         grad.beta3 = sum( dPredicted$treatment * invSigma.e),
                         grad.alpha1 = sum( bernoulliResidual ),
                         grad.alpha2 = sum( dPredicted$time * bernoulliResidual ),
                         grad.alpha3 = sum( dPredicted$treatment * bernoulliResidual ),
                         grad.gamma = sum( dPredicted$bDraw * bernoulliResidual ),
                         grad.sigma.b = grad.sigma.b,
                         grad.sigma = grad.sigma
                       )
                     }) %>% bind_rows %>% summarize_all(mean)
      
      out <- -2*unlist(temp, use.names = F)
      
      out
      
    }
    
    res <- optim(par=unlist(pars, use.names = F),
                 fn=minusTwoLogLikelihood,
                 gr=minusTwoScore,
                 method = 'L-BFGS-B',
                 control = list(factr=1e9),
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
    
    if(iter < 5)
    {
      mcSamples <- iter
    } else {
      mcSamples <- as.integer(min(mcSamples * 1.5, max(1e7 / nrow(d), 250))) #increase the sample size slowly
    }
    iter <- iter + 1
  }
  
  flog.trace(paste0('Sample ', key, ': EM result for spm: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ','), ' , deviance = ', format(currentMinus2LL, digits=7) ) )
  
  x0 <- unlist(pars)
  # hh <- hessian(function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5])), x0[c(1,2,3,8,9)])
  hh <- optimHess(par = x0[c(1,2,3,8,9)], 
                  fn = function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5])),
                  gr = function(x) minusTwoScore(c(x[1:3], x0[4:7], x[4:5]))[c(1,2,3,8,9)])
  
  out <- c(key, pars$beta, pars$sigma.b, pars$sigma, 2*diag(solve(hh)) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}
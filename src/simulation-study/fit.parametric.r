fit.parametric <- function(d, key) {
  flog.debug(paste0('Fitting parametric model to sample ', key))
  
  #needed?
  missingIds <- which(d$r == 0)
  
  #needed?
  # X <- model.matrix(as.formula(~treatment+time), d)
  
  #Fit ignorable model to find initial values for parameters
  
  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  pars <- list(beta = fixef(m),
               alpha = c(0,0,0),
               gamma = 0,
               sigma = sigma(m),
               sigma.b = as.data.frame(VarCorr(m))$sdcor[1])
  
  ##Loop for stochastic-EM
  previousMinus2LL <- Inf
  currentMinus2LL <- Inf
  iter <- 1
  minusTwoLogLikelihood <- NA
  
  while (coalesce(abs(previousMinus2LL-currentMinus2LL), Inf) > 1e-6 && iter <= 10000) {
    
    flog.trace(paste('SEM iteration', iter))
    
    # Stochastic expectation
    dPredicted <- d %>% mutate(
      bDraw = rep(rnorm(n=length(unique(d$subject)), sd=pars$sigma.b), each = 5),
      eDraw = rnorm(n=nrow(d), sd=pars$sigma),
      yPred = case_when(
        r == 1 ~ y,
        r == 0 ~ bDraw + pars$beta[1] + pars$beta[2] * time + pars$beta[3] * treatment + eDraw
      ))
    
    # Minimization step
    minusTwoLogLikelihood <- function(x) {
      
      if(length(x) != 9) flog.error('Incorrect input in logLikelihood for parametric model')
      
      beta <- x[1:3]
      alpha <- x[4:6]
      gamma <- x[7]
      sigma <- x[8]
      sigma.b <- x[9]
      
      logfyiri <- function(dd, key2) {
        yi <- dd$yPred
        ri <- dd$r
        Xi <- dd %>% mutate(intercept=1) %>% select(intercept, time, treatment) %>% as.matrix
        
        
        fyiri.conditional <- function(bi) {
          temp <- dmvnorm( yi, mean = Xi %*% beta + bi, sigma = sigma*diag(length(yi)) ) *
            prod(dbinom( ri, size = 1, prob = plogis(Xi %*% alpha + gamma * bi )))
          temp
        }
        
        #estimate subject log-likelihood contribution via Gauss-Hermite w/ Laplace approximation
        temp <- -2*data.frame(contrib=
                     log(aghQuad(g = function(bi) sapply(bi, fyiri.conditional), 
                                 muHat = 0, sigmaHat = sigma.b, 
                                 rule = gaussHermiteData(10))))
        
        if(!is.finite(temp$contrib[1])) {
          print(x)
        }
        
        temp
      }
      
      #sum the contributions
      dPredicted %>% group_by(subject) %>% group_modify(logfyiri) %>% pull(contrib) %>% sum
    }
    
    #Minimize the -2LL
    res <- optim(unlist(pars), minusTwoLogLikelihood, 
                 method = c('L-BFGS-B'),
                 lower = c(rep(-Inf, 7), 1e-4, 1e-4), #the 1e-4 are positivity constraints on variances
                 upper = Inf,
                 hessian = FALSE,
                 control = list(trace=1, REPORT=1, maxit=10)
    )
    if(res$convergence > 1) flog.error('Parametric model likelihood did not converge')
    
    pars <- list(beta = res$par[1:3],
                 alpha = res$par[4:6],
                 gamma = res$par[7],
                 sigma = res$par[8],
                 sigma.b = res$par[9])
    
    previousMinus2LL <- currentMinus2LL
    currentMinus2LL <- res$value
    iter <- iter + 1
  }
  
  flog.trace('Computing Hessian...')
  hh <- hessian(minusTwoLogLikelihood, unlist(pars))
  out <- c(pars$beta, pars$sigma.b, pars$sigma, diag(solve(hh))[c(1,2,3,9,8)] )
  names(out) <- c('intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}
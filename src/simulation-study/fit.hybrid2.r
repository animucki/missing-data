#Fit the class model of Lin, McCulloch, and Rosenheck, but with Bernoulli missingness
fit.hybrid2 <- function(d) {
  key <- as.integer(d$sample[1])
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
  
  integrationRule <- gaussHermiteData(20) # numerical rule for the Gauss-Hermite integral
  
  flog.trace(paste0('Sample ', key, ': direct hybrid initial pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ','),
                    ', deviance = ', format(res$value, digits=7)) )
  
  # Define likelihood function
  minusTwoLogLikelihood <- function(x) {
    
    beta <- x[1:3]
    sigma.b <- x[4]
    sigma <- x[5]
    alpha <- x[6:8]
    gamma <- x[9]
    delta <- x[10]
    mu <- x[10 + 1:(nClasses-1)]
    eta <- x[10 + (nClasses-1) + 1:(nClasses-1)]
    
    ll <- d %>% group_split(subject) %>% lapply(function(di) {
      
      ni <- sum(!is.na(di$y)) #convenience: number of observations for subject i
      cik <- softmax(c(0,eta))
      
      ghQuad(
        f=function(bi) {
          unlist(lapply(bi, function(bi) {
            likelihoodOverClasses <- 0
            classProb <- softmax(eta)
            
            for(k in 1:length(classProb)) {
              likelihoodOverClasses <- likelihoodOverClasses +
                prod(dnorm(
                  x= di$y - bi - beta[1] - di$time * beta[2] - di$treatment * beta[3] - mu[k],
                  sd = sigma
                ), na.rm = TRUE) *
                prod(dbinom(
                  x = di$r,
                  size = 1,
                  prob = plogis(alpha[1] + di$time * alpha[2] + di$treatment * alpha[3] + gamma * bi + delta * mu[k])
                )) *
                classProb[k]
            }
            return(likelihoodOverClasses)
          })) *
            dnorm(bi, sd=sigma.b) / exp(-bi^2) #the exponential factor is because fastGHquad multiplies f by exp(-bi^2)
        },
        rule=integrationRule
      ) %>% log
      
    }) %>% unlist %>% sum
    
    -2*ll
  }
  
  res <- optim(par=unlist(pars, use.names = F), 
               fn=minusTwoLogLikelihood,
               method = 'L-BFGS-B',
               control = list(trace=1, REPORT=1),
               lower = c(rep(-Inf, 7), 1e-4, 1e-4, rep(-Inf, 2*(nClasses-1))),
               upper = Inf,
               hessian = FALSE
  )
  
  if(res$convergence > 0) flog.error(paste('hybrid likelihood did not converge, code',res$convergence))
  
  pars <- list(beta = res$par[1:3],
               sigma.b = res$par[4],
               sigma = res$par[5],
               alpha = res$par[6:8],
               gamma = res$par[9],
               delta = res$par[10],
               mu = res$par[10 + 1:(nClasses-1)],
               eta = res$par[10 + (nClasses-1) + 1:(nClasses-1)])

  flog.trace(paste0('Sample ', key, ': EM result for hybrid: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ','), ' , deviance = ', format(res$value, digits=7) ) )
  
  x0 <- unlist(pars)
  hh <- hessian(function(x) minusTwoLogLikelihood(c(x[1:5], x0[-(1:5)])), x0[1:5])
  
  print(hh)
  
  out <- c(key, pars$beta, pars$sigma.b, pars$sigma, 2*diag(solve(hh)) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  
  flog.trace(paste0('Sample ', key, ': hybrid complete, incl. Hessian.'))
  
  as.data.frame(as.list(out))
}
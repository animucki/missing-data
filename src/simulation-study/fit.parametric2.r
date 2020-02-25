#Fit the shared-parameter model, but deterministically, using numerical integration over the random effects

fit.parametric2 <- function(d) {
  
  key <- as.integer(d$sample[1])
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
  
  integrationRule <- gaussHermiteData(20) # numerical rule for the Gauss-Hermite integral
  
  flog.trace(paste0('Sample ', key, ': direct algorithm initial spm pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ',')))
  
  # Formulate likelihood function
  minusTwoLogLikelihood <- function(x) {
    beta <- x[1:3]
    alpha <- x[4:6]
    gamma <- x[7]
    sigma.b <- x[8]
    sigma <- x[9]
    
    ll <- d %>% group_split(subject) %>% lapply(function(di) {
      
      ni <- sum(!is.na(di$y)) #convenience: number of observations for subject i
      
      ghQuad(
        f=function(bi) {
          unlist(lapply(bi, function(bi) {
            prod(dnorm(
              x= di$y - bi - beta[1] - di$time * beta[2] - di$treatment * beta[3] ,
              sd = sigma
            ), na.rm = TRUE) *
              prod(dbinom(
                x = di$r,
                size = 1,
                prob = plogis(alpha[1] + di$time * alpha[2] + di$treatment * alpha[3] + gamma * bi)
              ))
          })) *
            dnorm(bi, sd=sigma.b) / exp(-bi^2) #the exponential factor is because fastGHquad multiplies f by exp(-bi^2)
        },
        rule=integrationRule
      ) %>% log
      
    }) %>% unlist %>% sum
    
    -2*ll
  }
  
  #Run EM
  res <- optim(par=unlist(pars, use.names = F),
               fn=minusTwoLogLikelihood,
               method = 'L-BFGS-B',
               control = list(trace = 1, maxit = 2000),
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
  
  flog.trace(paste0('Sample ', key, ': direct result for spm: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ','),
                    ', deviance = ', format(res$value, digits=7)) )
  
  x0 <- unlist(pars)
  hh <- hessian(function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5])), x0[c(1,2,3,8,9)])
  
  out <- c(key, pars$beta, pars$sigma.b, pars$sigma, 2*diag(solve(hh)) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  
  flog.trace(paste0('Sample ', key, ': spm complete, incl. Hessian.'))
  
  as.data.frame(as.list(out))
}
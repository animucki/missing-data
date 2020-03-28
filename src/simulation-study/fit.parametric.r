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
  
  ##Loop for MCEM
  previousPars <- Inf
  currentPars <- unlist(pars)
  
  iter <- 1
  
  mcSamples <- 5
  
  minusTwoLogLikelihood <- NA

  dPredictedList <- list()

  nTimesCriterionMet <- 0
  crit <- Inf

  while (nTimesCriterionMet < 3 && iter <= 100) {
    
    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ', MC samples: ', mcSamples,', pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ','),
                      ', crit = ', format(crit, digits=7) ))

    # Monte Carlo Expectation step
    # Initialize the list of simulated datasets
    for (mc in 1:mcSamples) {
      dPredictedList[[mc]]  <- d
      dPredictedList[[mc]]$bDraw <- NA_real_
    }

    for (i in seq_along(unique(d$subject))) {
      # simulate

      di <- d %>% filter(subject==i)

      bVec <- arms(n_samples = mcSamples, lower = -10, upper = 10, metropolis = FALSE,
                   log_pdf = function(bi) {
                     sum(dnorm(
                       x=di$y,
                       mean=pars$beta[1] + di$time * pars$beta[2] + di$treatment * pars$beta[3] + bi,
                       sd = pars$sigma,
                       log = T), na.rm = T)+ #na.rm=T skips the missing observations
                       sum(dbinom(
                         x=di$r,
                         size = 1,
                         prob = plogis(pars$alpha[1] + di$time * pars$alpha[2] + di$treatment * pars$alpha[3] + pars$gamma * bi ),
                         log = T)[-1]) + #the [-1] skips the baseline observation
                       dnorm(
                         x=bi,
                         sd = pars$sigma.b,
                         log = T)
                   })

      # save simulated values
      for (mc in 1:mcSamples) {
        dPredictedList[[mc]]$bDraw[d$subject==i] <- bVec[mc]
      }
    }

    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ', E step completed'))

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
                         sum(dnorm(
                           x=dPredicted$y,
                           mean=beta[1] + Xtemp %*% beta[c(2,3)] + dPredicted$bDraw, 
                           sd = sigma, 
                           log = T), na.rm = T)+ #na.rm=T skips the missing observations
                           sum(dbinom(
                             x=dPredicted$r, 
                             size = 1, 
                             prob = plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] + gamma * dPredicted$bDraw ),
                             log = T)[-seq(from = 1, to = nrow(d), by=nTimePoints)]) + #the seq(...) expression skips the baseline observations
                           sum(dnorm(
                             x=dPredicted$bDraw[seq(1, nrow(dPredicted), by = nTimePoints)], 
                             sd = sigma.b, 
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
                       normalResidual <- dPredicted$y - beta[1] - Xtemp %*% beta[c(2,3)] - dPredicted$bDraw
                       bernoulliResidual <- dPredicted$r - plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] + gamma * dPredicted$bDraw )
                       data.frame(
                         grad.beta1 = sum( normalResidual, na.rm = TRUE )/sigma^2,
                         grad.beta2 = sum( dPredicted$time *  normalResidual, na.rm = TRUE )/sigma^2,
                         grad.beta3 = sum( dPredicted$treatment * normalResidual, na.rm = TRUE )/sigma^2,
                         grad.alpha1 = sum( (bernoulliResidual)[-seq(1, nrow(dPredicted), by = nTimePoints)] ),
                         grad.alpha2 = sum( (dPredicted$time * bernoulliResidual)[-seq(1, nrow(dPredicted), by = nTimePoints)] ),
                         grad.alpha3 = sum( (dPredicted$treatment * bernoulliResidual)[-seq(1, nrow(dPredicted), by = nTimePoints)] ),
                         grad.gamma = sum( (dPredicted$bDraw * bernoulliResidual)[-seq(1, nrow(dPredicted), by = nTimePoints)] ),
                         grad.sigma.b = sum( (dPredicted$bDraw[seq(1, nrow(dPredicted), by = nTimePoints)]^2 - sigma.b^2 )/sigma.b^3),
                         grad.sigma = sum( normalResidual^2 - sigma^2, na.rm = TRUE )/sigma^3)
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
  
  flog.trace(paste0('Sample ', key, ': EM result for spm: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ',')) )
  
  x0 <- unlist(pars)
  # hh <- hessian(function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5])), x0[c(1,2,3,8,9)])
  hh <- optimHess(par = x0[c(1,2,3,8,9)], 
                  fn = function(x) minusTwoLogLikelihood(c(x[1:3], x0[4:7], x[4:5])),
                  gr = function(x) minusTwoScore(c(x[1:3], x0[4:7], x[4:5]))[c(1,2,3,8,9)])
  
  out <- c(key, pars$beta, pars$sigma.b, pars$sigma, sqrt(diag(2*solve(hh))) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}
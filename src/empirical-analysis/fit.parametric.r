#Fit the shared-parameter model
fit.parametric <- function(y, r, X, W, init) {

  ni <- ncol(y)
  n <- nrow(y)
  p <- ncol(X)
  pr <- ncol(W)

  pars <- list(beta = init[1:p],
               alpha = rep(0,pr),
               gamma = 0,
               sigma.b = init[p+1],
               sigma = init[p+2]
               )

  ##Loop for MCEM
  currentPars <- unlist(pars)

  iter <- 1

  mcSamples <- 10

  minusTwoLogLikelihood <- NA

  dPredictedList <- list()

  nTimesCriterionMet <- 0
  crit <- Inf

  while (nTimesCriterionMet < 3 && iter <= 100) {

    flog.trace(paste0('EM iteration ', iter, ', MC samples: ', mcSamples,', pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ','),
                      ', crit = ', format(crit, digits=7) ))

    # Monte Carlo Expectation step
    # Initialize the list of simulated random intercepts
    for (mc in 1:mcSamples) {
      dPredictedList[[mc]]  <- list()
      dPredictedList[[mc]]$bDraw <- NA_real_
    }

    for (i in 1:n) {
      # simulate

      bVec <- arms(n_samples = mcSamples, lower = -10, upper = 10, metropolis = FALSE,
                   log_pdf = function(bi) {
                     sum(dnorm(
                       x = y[i,],
                       mean = X[5*(i-1) + 1:5,] %*% pars$beta + bi,
                       sd = pars$sigma,
                       log = T), na.rm = T)+ #na.rm=T skips the missing observations
                       sum(dbinom(
                         x=r[i,],
                         size = 1,
                         prob = plogis( W[5*(i-1) + 1:5,] %*% pars$alpha + pars$gamma * bi ),
                         log = T)) +
                       dnorm(
                         x=bi,
                         sd = pars$sigma.b,
                         log = T)
                   })

      # save simulated values
      for (mc in 1:mcSamples) {
        dPredictedList[[mc]]$bDraw[i] <- bVec[mc]
      }
    }

    flog.trace(paste0('EM iteration ', iter, ', E step completed'))

    # Minimization step
    minusTwoLogLikelihood <- function(x) {

      beta <- x[1:p]
      alpha <- x[p + 1:pr]
      gamma <- x[p + pr + 1]
      lsigma.b <- x[p + pr + 2]
      sigma.b <- exp(lsigma.b)
      lsigma <- x[p + pr + 3]
      sigma <- exp(lsigma)

      temp <- lapply(dPredictedList,
                       function(dPredicted) {
                         sum(dnorm(
                           x = as.vector(t(y)),
                           mean= X %*% beta + rep(dPredicted$bDraw, each = ni),
                           sd = sigma,
                           log = T), na.rm = T)+ #na.rm=T skips the missing observations
                           sum(dbinom(
                             x= as.vector(t(r)),
                             size = 1,
                             prob = plogis( W %*% alpha + gamma * rep(dPredicted$bDraw, each = ni) ),
                             log = T)) +
                           sum(dnorm(
                             x=dPredicted$bDraw,
                             sd = sigma.b,
                             log = T))
                       })

      out <- -2*mean(unlist(temp))

      out

    }

    minusTwoScore <- function(x) {

      beta <- x[1:p]
      alpha <- x[p + 1:pr]
      gamma <- x[p + pr + 1]
      lsigma.b <- x[p + pr + 2]
      sigma.b <- exp(lsigma.b)
      lsigma <- x[p + pr + 3]
      sigma <- exp(lsigma)

      temp <- lapply(dPredictedList,
                     function(dPredicted) {
                       normalResidual <- as.vector(t(y)) - X %*% beta - rep(dPredicted$bDraw, each = ni)
                       bernoulliResidual <- as.vector(t(r)) - plogis(W %*% alpha + gamma * rep(dPredicted$bDraw, each = ni) )

                       grad <- data.frame(0)
                       for(j in 1:p) {
                         grad[[paste0('grad.beta',j)]] <- sum( X[,j] * normalResidual, na.rm = TRUE )/sigma^2
                       }
                       for(j in 1:pr) {
                         grad[[paste0('grad.alpha',j)]] <- sum( W[,j] * bernoulliResidual )
                       }

                       grad$grad.gamma = sum( (rep(dPredicted$bDraw, each = ni) * bernoulliResidual))
                       grad$grad.sigma.b = sum( (dPredicted$bDraw^2 - sigma.b^2 )/sigma.b^2)
                       grad$grad.sigma = sum( normalResidual^2 - sigma^2, na.rm = TRUE )/sigma^2

                       return(grad[,-1])
                     }) %>% bind_rows %>% summarize_all(mean)

      out <- -2*unlist(temp, use.names = F)

      out

    }

    res <- optim(par=unlist(pars, use.names = F),
                 fn=minusTwoLogLikelihood,
                 gr=minusTwoScore,
                 # method = 'L-BFGS-B',
                 method = 'BFGS',
                 # lower = c(rep(-Inf, p+pr+1), 1e-4, 1e-4),
                 # upper = Inf,
                 hessian = FALSE
                 )

    if(res$convergence > 1) flog.error(paste('Parametric model likelihood did not converge, code',res$convergence))

    pars <- list(beta = res$par[1:p],
                 alpha = res$par[p + 1:pr],
                 gamma = res$par[p + pr + 1],
                 sigma.b = exp(res$par[p + pr + 2]),
                 sigma = exp(res$par[p + pr + 3])
    )

    previousPars <- currentPars
    currentPars <- unlist(pars)


    mcSamples <- as.integer(min(mcSamples * 1.3, 500)) #increase the sample size slowly
    iter <- iter + 1

    #stopping criteria calculation
    crit <- coalesce(mean( (previousPars-currentPars)^2 )/(mean(previousPars^2)+1e-3), Inf)
    if(crit<1e-4) {
      nTimesCriterionMet <- nTimesCriterionMet + 1
    } else {
      nTimesCriterionMet <- 0
    }

  }

  flog.trace(paste0('EM result for spm: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ',')) )

  hess <- hessian(func = minusTwoLogLikelihood,
                  x = unlist(pars, use.names = F))

  return(list(res = res$value, pars = pars, hess = hess))
}
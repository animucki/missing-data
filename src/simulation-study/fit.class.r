#Fit the class model of Lin, McCulloch, and Rosenheck, but with Bernoulli missingness
fit.class <- function(d) {
  key <- as.integer(d$sample[1])
  set.seed(4114L + key)
  flog.debug(paste0('Fitting class model to sample ', key))

  nSubjects <- length(unique(d$subject))
  nTimePoints <- d %>% group_by(subject) %>% summarize(n=n(), .groups="drop") %>% pull(n) %>% max
  nClasses <- 3

  # prepare grad0 (initial gradient)
  grad0 <- data.frame(grad.beta1 = 0, grad.beta2 = 0, grad.beta3 = 0,
                      grad.alpha1 = 0, grad.alpha2 = 0, grad.alpha3 = 0,
                      grad.sigma.b = 0, grad.sigma = 0, grad.delta = 0)
  for (k in 2:nClasses) {
    grad0[[paste0('grad.mu',k)]] <- 0
  }
  for (k in 2:nClasses) {
    grad0[[paste0('grad.eta',k)]] <- 0
  }
  
  #Fit ignorable model to find initial values for parameters
  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  pars <- list(beta = fixef(m),
               alpha = c(0,0,0),
               lsigma.b = log(as.data.frame(VarCorr(m))$sdcor[1]),
               lsigma = log(sigma(m)),
               delta = 0,
               mu = rep(0, nClasses-1),
               eta = rep(0, nClasses-1)
  )

  # smart initialization to make separation of class intercepts more likely
  pars$mu <- c(-1,1) * exp(pars$lsigma.b)

  currentPars <- unlist(pars)

  iter <- 1

  mcSamples <- 5

  minusTwoLogLikelihood <- NA

  dPredictedList <- list()

  nTimesCriterionMet <- 0
  crit <- Inf

  dList <- d %>% group_split(subject)
  doList <- d %>% filter(r==1) %>% group_split(subject)

  while (nTimesCriterionMet < 3 && iter <= 100) {

    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ', MC samples: ', mcSamples,', pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ','),
                      ', crit = ', format(crit, digits = 4)) )

    # Monte Carlo Expectation step
    # Initialize the list of simulated datasets
    for (mc in 1:mcSamples) {
      dPredictedList[[mc]] <- matrix(NA_integer_, nrow = nSubjects, ncol = nClasses)
    }

    #simulate
    for (i in 1:nSubjects) {

      # find conditional distribution of ci for this subject
      pr <- rep(NA_real_, nClasses)
      for (k in 1:nClasses) {
        cik <- rep(0, nClasses)
        cik[k] <- 1

        pr[k] <- exp(dmultinormal(x = doList[[i]]$y,
                                  mean = pars$beta[1] + pars$beta[2] * doList[[i]]$time + pars$beta[3] * doList[[i]]$treatment + as.vector(c(0,pars$mu) %*% cik),
                                  sigma = as.vector( exp(pars$lsigma)^2 * diag(nrow(doList[[i]])) + exp(pars$lsigma.b)^2),
                                  log = TRUE) +
                       sum( dbinom(x = dList[[i]]$r,
                                   prob = plogis(pars$alpha[1] + dList[[i]]$time * pars$alpha[2] + dList[[i]]$treatment * pars$alpha[3] +
                                                   pars$delta * c(0, pars$mu)[k]),
                                   size = 1,
                                   log = TRUE)[-1] ) +
                       dmultinomial( x = cik,
                                     size = 1,
                                     prob = softmax(c(0,pars$eta)),
                                     log = TRUE))
      }

      ci <- rmultinomial(n = mcSamples, size = 1, prob = pr)

      # Save values
      for (mc in 1:mcSamples) {
        dPredictedList[[mc]][i,] <- ci[mc,]
      }
    }

    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ' E step completed'))

    # Minimization step
    minusTwoLogLikelihood <- function(x) {

      beta <- x[1:3]
      alpha <- x[4:6]
      lsigma.b <- x[7]
      lsigma <- x[8]
      delta <- x[9]
      mu <- x[9 + 1:(nClasses-1)]
      eta <- x[9 + (nClasses-1) + 1:(nClasses-1)]

      sigma.b <- exp(lsigma.b)
      sigma <- exp(lsigma)

      tmp <- lapply(dPredictedList,
                    function(dObj) {
                      ll <- 0
                      for (i in 1:nSubjects) {
                        ll <- ll +
                          dmultinormal(x = doList[[i]]$y,
                                       mean = beta[1] + beta[2] * doList[[i]]$time + beta[3] * doList[[i]]$treatment + as.vector(c(0,mu) %*% dObj[i,]),
                                       sigma = as.vector( sigma^2 * diag(nrow(doList[[i]])) + sigma.b^2 ),
                                       log = T) +
                          sum( dbinom(x = dList[[i]]$r,
                                      prob = plogis(alpha[1] + dList[[i]]$time * alpha[2] + dList[[i]]$treatment * alpha[3] +
                                                      delta * as.vector(c(0,mu) %*% dObj[i,]) ),
                                      size = 1,
                                      log = TRUE)[-1] )
                      }
                      ll <- ll +
                        sum(dmultinomial( x = dObj,
                                          size = 1,
                                          prob = softmax(c(0,eta)),
                                          log = T))

                      return(ll)

                    })

      out <- -2*mean(unlist(tmp))

      out

    }

    minusTwoScore <- function(x) {

      beta <- x[1:3]
      alpha <- x[4:6]
      lsigma.b <- x[7]
      lsigma <- x[8]
      delta <- x[9]
      mu <- x[9 + 1:(nClasses-1)]
      eta <- x[9 + (nClasses-1) + 1:(nClasses-1)]

      sigma.b <- exp(lsigma.b)
      sigma <- exp(lsigma)

      out <- lapply(dPredictedList,
                    function(dObj) {

                      grad <- grad0
                      for (i in 1:nSubjects) {
                        #precalculate quantities
                        ni <- nrow(doList[[i]])
                        Xi <- cbind(
                          matrix(c(rep(1,ni),
                                   doList[[i]]$time,
                                   doList[[i]]$treatment),
                                 nrow = ni),
                          matrix(rep(dObj[i,-1],
                                     each = ni),
                                 nrow = ni))

                        invSigma <- solve( sigma^2 * diag(ni) + sigma.b^2 )
                        ei <- doList[[i]]$y - Xi %*% c(beta, mu)

                        tau <- invSigma %*% (ei %*% t(ei) - (sigma^2 * diag(ni) + sigma.b^2)) %*% invSigma

                        # grad.beta and grad.mu can be calculated at once (score equations as derived by Jennrich and Schluchter)
                        grad[,c(1,2,3,9 + 1:(nClasses-1))] <- grad[,c(1,2,3,9 + 1:(nClasses-1))] + t(Xi) %*% invSigma %*% ei
                        # then, we have the variance components

                        grad[,7] <- grad[,7] + 1/2 * tr(tau %*% (matrix(2*sigma.b, nrow = ni, ncol = ni)))
                        grad[,8] <- grad[,8] + 1/2 * tr(tau %*% (2*sigma * diag(ni)))

                      }

                      grad[,7] <- grad[,7] * sigma.b
                      grad[,8] <- grad[,8] * sigma

                      #grad.alpha
                      bernoulliResidual <- d$r - plogis(alpha[1] + d$time * alpha[2] + d$treatment * alpha[3] +
                                                          delta * as.vector(dObj %*% c(0, mu)))
                      grad[4] <- sum(bernoulliResidual[-seq(from = 1, to = nrow(d), by = nTimePoints)])
                      grad[5] <- sum((d$time * bernoulliResidual)[-seq(from = 1, to = nrow(d), by = nTimePoints)])
                      grad[6] <- sum((d$treatment * bernoulliResidual)[-seq(from = 1, to = nrow(d), by = nTimePoints)])

                      #grad.delta
                      grad[9] <- sum((rep(as.vector(dObj %*% c(0, mu)), each = nTimePoints) *
                        bernoulliResidual)[-seq(from = 1, to = nrow(d), by = nTimePoints)])

                      #grad.eta
                      grad[,9 + (nClasses-1) + 1:(nClasses-1)] <- colSums(dObj[,-1] - rep(softmax(c(0,eta))[-1], each = nrow(dObj)))

                      return(grad)
                    }) %>% bind_rows %>% summarize_all(mean) %>% unlist

      -2*out
    }

    #lowerBounds <- rep(-Inf, length(unlist(pars)))
    #lowerBounds[c(6,7,8, 8 + 2*(nClasses-1) + 1:nClasses, length(unlist(pars)))] <- 1e-4 #the three variances, and the lambdas, have to be positive

    #ctrl <- list(trace = 6, REPORT = 1)
    ctrl <- list(maxit=min(25*iter,100), reltol = 1e-4)
    if(nTimesCriterionMet == 2) ctrl <- list(maxit=250, reltol = 1e-8)

    res <- optim(par=unlist(pars, use.names = F),
                 fn=minusTwoLogLikelihood,
                 gr=minusTwoScore,
                 method = 'BFGS',
                 control = ctrl,
                 #lower = lowerBounds,
                 #upper = Inf,
                 hessian = FALSE
    )

    if(res$convergence == 1 && nTimesCriterionMet >= 2) flog.warn('class likelihood did not converge within 250 iterations in final step')
    if(res$convergence > 0) flog.warn(paste('Class model likelihood did not converge, code',res$convergence))

    pars <- list(beta = res$par[1:3],
                 alpha = res$par[4:6],
                 lsigma.b = res$par[7],
                 lsigma = res$par[8],
                 delta = res$par[9],
                 mu = res$par[9 + 1:(nClasses-1)],
                 eta = res$par[9 + (nClasses-1) + 1:(nClasses-1)]
    )

    previousPars <- currentPars
    currentPars <- unlist(pars)

    mcSamples <- floor(min(mcSamples + 5, 250)) #increase the sample size slowly
    iter <- iter + 1

    #stopping criteria calculation
    crit <- coalesce(mean( (previousPars-currentPars)^2 )/(mean(previousPars^2)+1e-3), Inf)
    if(crit<5e-4) {
      nTimesCriterionMet <- nTimesCriterionMet + 1
    } else {
      nTimesCriterionMet <- 0
    }
  }

  flog.trace(paste0('Sample ', key, ': EM result for class: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ',')))

  x0 <- unlist(pars)
  #hh1 <- optimHess(par = x0,
  #                fn = minusTwoLogLikelihood,
  #                gr = minusTwoScore
  #)

  hh1 <- numDeriv::jacobian(func = minusTwoScore,
                            x = x0)

  ainv <- rep(1,length(x0))
  ainv[7] <- exp(pars$lsigma.b)
  ainv[8] <- exp(pars$lsigma)
  gradi <- rep(0, length(x0))
  grads <- minusTwoScore(x0)
  gradi[7] <- grads[7]/exp(pars$lsigma.b)^2
  gradi[8] <- grads[8]/exp(pars$lsigma)^2

  hh <- diag(ainv) %*% hh1 %*% diag(ainv) + diag(gradi)

  out <- c(key, pars$beta, exp(pars$lsigma.b), exp(pars$lsigma), sqrt(diag(2*solve(hh))[c(1,2,3,7,8)]) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}
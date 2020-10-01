#Fit the Tseng et al. model
fit.tseng <- function(d) {
  #IMPORTANT: this function assumes that the observations are ordered by time!!!

  key <- as.integer(d$sample[1])
  set.seed(88810L + key)
  flog.debug(paste0('Fitting parametric model to sample ', key))

  nTimePoints <- d %>% group_by(subject) %>% summarize(n=n()) %>% pull(n) %>% max

  #Fit ignorable model to find initial values for parameters

  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  pars <- list(beta = fixef(m),
               alpha = c(0,0,0),
               lsigma.b1 = log(as.data.frame(VarCorr(m))$sdcor[1]),
               lsigma.b2 = log(as.data.frame(VarCorr(m))$sdcor[1]),
               tau.b = 0,
               lsigma = log(sigma(m))
  )

  ##Loop for MCEM
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
      dPredictedList[[mc]]$b1Draw <- NA_real_
      dPredictedList[[mc]]$b2Draw <- NA_real_
    }

    for (i in seq_along(unique(d$subject))) {
      # simulate

      di <- d %>% filter(subject==i)

      bMat <- dlm::arms(n.sample = mcSamples,
                        y.start = c(0,0),
                        indFunc = function(bi) { sum(bi^2) < 10 },
                        myldens = function(bi) {
                          sum(dnorm(
                            x=di$y,
                            mean=pars$beta[1] + di$time * pars$beta[2] + di$treatment * pars$beta[3] + bi[1],
                            sd = exp(pars$lsigma),
                            log = T), na.rm = T)+ #na.rm=T skips the missing observations
                            sum(dbinom(
                              x=di$r,
                              size = 1,
                              prob = plogis(pars$alpha[1] + di$time * pars$alpha[2] + di$treatment * pars$alpha[3] + pars$gamma * bi[2] ),
                              log = T)[-1]) + #the [-1] skips the baseline observation
                            dmultinormal(
                              x = bi,
                              sigma = c(exp(2*pars$lsigma.b1), pars$tau.b, pars$tau.b, exp(2*pars$lsigma.b2)),
                              log = T)
                        }
      )

      # save simulated values
      for (mc in 1:mcSamples) {
        dPredictedList[[mc]]$b1Draw[d$subject==i] <- bMat[mc,1]
        dPredictedList[[mc]]$b2Draw[d$subject==i] <- bMat[mc,2]
      }

    }

    flog.trace(paste0('Sample ', key, ': EM iteration ', iter, ', E step completed'))

    Xtemp <- as.matrix(d %>% select(time, treatment))
    # Minimization step
    minusTwoLogLikelihood <- function(x) {

      if(length(x) != 10) flog.error('Incorrect input in logLikelihood for parametric model')

      beta <- x[1:3]
      alpha <- x[4:6]
      lsigma.b1 <- x[7]
      lsigma.b2 <- x[8]
      tau.b <- x[9]
      lsigma <- x[10]

      sigma.b1 <- exp(lsigma.b1)
      sigma.b2 <- exp(lsigma.b2)
      sigma <- exp(lsigma)

      temp <- lapply(dPredictedList,
                     function(dPredicted) {
                       sum(dnorm(
                         x=dPredicted$y,
                         mean=beta[1] + Xtemp %*% beta[-1] + dPredicted$b1Draw,
                         sd = sigma,
                         log = T), na.rm = T)+ #na.rm=T skips the missing observations
                         sum(dbinom(
                           x=dPredicted$r,
                           size = 1,
                           prob = plogis(alpha[1] + Xtemp %*% alpha[-1] + dPredicted$b2Draw ),
                           log = T)[-seq(from = 1, to = nrow(d), by=nTimePoints)]) + #the seq(...) expression skips the baseline observations
                         sum(dmultinormal(
                           x=matrix(c(dPredicted$b1Draw[seq(1, nrow(dPredicted), by = nTimePoints)],
                                      dPredicted$b2Draw[seq(1, nrow(dPredicted), by = nTimePoints)]),
                                    ncol = 2),
                           sigma = c(sigma.b1^2, tau.b, tau.b, sigma.b2^2),
                           log = T))
                     })

      out <- -2*mean(unlist(temp))

      out

    }

    minusTwoScore <- function(x) {

      if(length(x) != 10) flog.error('Incorrect input in score for parametric model')

      beta <- x[1:3]
      alpha <- x[4:6]
      lsigma.b1 <- x[7]
      lsigma.b2 <- x[8]
      tau.b <- x[9]
      lsigma <- x[10]

      sigma.b1 <- exp(lsigma.b1)
      sigma.b2 <- exp(lsigma.b2)
      sigma <- exp(lsigma)

      sigma_invt <- t(solve(matrix(c(sigma.b1^2, tau.b, tau.b, sigma.b2^2), ncol = 2)))

      temp <- lapply(dPredictedList,
                     function(dPredicted) {
                       normalResidual <- dPredicted$y - beta[1] - Xtemp %*% beta[c(2,3)] - dPredicted$b1Draw
                       bernoulliResidual <- dPredicted$r - plogis(alpha[1] + Xtemp %*% alpha[c(2,3)] + dPredicted$b2Draw )
                       tempmc <- data.frame(
                         grad.beta1 = sum( normalResidual, na.rm = TRUE )/sigma^2,
                         grad.beta2 = sum( dPredicted$time *  normalResidual, na.rm = TRUE )/sigma^2,
                         grad.beta3 = sum( dPredicted$treatment * normalResidual, na.rm = TRUE )/sigma^2,
                         grad.alpha1 = sum( (bernoulliResidual)[-seq(1, nrow(dPredicted), by = nTimePoints)] ),
                         grad.alpha2 = sum( (dPredicted$time * bernoulliResidual)[-seq(1, nrow(dPredicted), by = nTimePoints)] ),
                         grad.alpha3 = sum( (dPredicted$treatment * bernoulliResidual)[-seq(1, nrow(dPredicted), by = nTimePoints)] ),
                         grad.sigma.b1 = 0,
                         grad.sigma.b2 = 0,
                         grad.tau.b = 0,
                         grad.sigma = sum( normalResidual^2 - sigma^2, na.rm = TRUE )/sigma^2)

                       for(i in seq(from = 1, to=nrow(dPredicted), by = nTimePoints)) {
                         u <- c(dPredicted$b1Draw[i], dPredicted$b2Draw[i])
                         tempmc[,c('grad.sigma.b1','grad.sigma.b2','grad.tau.b')] <-
                           tempmc[,c('grad.sigma.b1','grad.sigma.b2','grad.tau.b')] -
                             (sigma_invt %*% (u %o% u) %*% sigma_invt)[c(1,4,2)] -
                             sigma_invt[c(1,4,2)]
                       }

                       tempmc[,c('grad.sigma.b1','grad.sigma.b2')] <-
                         tempmc[,c('grad.sigma.b1','grad.sigma.b2')] *c(2,2)

                       tempmc
                     }) %>% bind_rows %>% summarize_all(mean)

      out <- -2*unlist(temp, use.names = F)

      gradTrue <- grad(minusTwoLogLikelihood, x)
      print(out)
      print(gradTrue)
      print(out/gradTrue)

      out

    }

    res <- optim(par=unlist(pars, use.names = F),
                 fn=minusTwoLogLikelihood,
                 #gr=minusTwoScore,
                 method = 'BFGS',
                 control = list(trace=1, REPORT=1),
                 #lower = c(rep(-Inf, 7), 1e-4, 1e-4),
                 #upper = Inf,
                 hessian = FALSE
    )

    if(res$convergence > 0) {
      flog.warn(paste('Parametric model likelihood did not converge, code',res$convergence))
      out <- c(key, rep(NA_real_, 2 * (length(pars$beta) + 2)) )
      names(out) <- c('sample','intercept', 'time', 'treatment', 'lsigma.b', 'lsigma',
                      'se.intercept', 'se.time', 'se.treatment', 'se.lsigma.b', 'se.lsigma')
      return(as.data.frame(as.list(out)))
    }

    pars <- list(beta = res$par[1:3],
                 alpha = res$par[4:6],
                 lsigma.b1 = res$par[7],
                 lsigma.b2 = res$par[8],
                 tau.b = res$par[9],
                 lsigma = res$par[10]
    )

    previousPars <- currentPars
    currentPars <- unlist(pars)


    mcSamples <- as.integer(min(mcSamples + 5, 250)) #increase the sample size slowly
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
  hh1 <- optimHess(par = x0,
                   fn = minusTwoLogLikelihood,
                   gr = minusTwoScore)
  ainv <- rep(1,length(x0))
  ainv[length(ainv)-2] <- exp(pars$lsigma.b)
  ainv[length(ainv)-1] <- exp(pars$lsigma)
  gradi <- rep(0, length(x0))
  grads <- minusTwoScore(x0)
  gradi[length(ainv)-2] <- grads[length(ainv)-2]/exp(pars$lsigma.b)^2
  gradi[length(ainv)-1] <- grads[length(ainv)-1]/exp(pars$lsigma)^2

  hh <- diag(ainv) %*% hh1 %*% diag(ainv) + diag(gradi)

  out <- c(key, pars$beta, exp(pars$lsigma.b), exp(pars$lsigma), sqrt(diag(2*solve(hh))[c(1,2,3,8,9)]) )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}
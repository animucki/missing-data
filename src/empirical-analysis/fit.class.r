#Fit the class model of Lin, McCulloch, and Rosenheck
fit.class <- function(y, r, X, W, nClasses, init) {

  ni <- ncol(y)
  n <- nrow(y)
  p <- ncol(X)
  pr <- ncol(W)

  # prepare grad0 (initial gradient)
  grad0 <- data.frame(0)
  for (j in 1:p) {
    grad0[[paste0('grad.beta',j)]] <- 0
  }
  grad0[['grad.alpha_time']] <- 0
  for (j in 2:pr) {
    grad0[[paste0('grad.alpha',j)]] <- 0
  }

  grad0$grad.sigma.b <- 0
  grad0$grad.sigma <- 0
  grad0$grad.theta <- 0

  for (k in 2:nClasses) {
    grad0[[paste0('grad.mu',k)]] <- 0
  }
  for (k in 2:nClasses) {
    grad0[[paste0('grad.eta',k)]] <- 0
  }
  for (k in 1:nClasses) {
    grad0[[paste0('grad.lambda',k)]] <- 0
  }

  if(!("time" %in% colnames(W))) stop("Time must be included in the W matrix")
  time <- W[,"time"]
  W <- W[,-which(colnames(W)=="time")]

  pars <- list(beta = init[1:p],
               alpha_time = 1e-2,
               alpha = rep(0,pr-1), # time parameter of the hazard function - the likelihood has a removable discontinuity at alpha_time==0 so don't start too close to it
               sigma.b = init[p+1],
               sigma = init[p+2],
               theta = 1e-2, # variance of the frailty distribution
               mu = rep(0, nClasses-1),
               eta = rep(0, nClasses-1),
               lambda = rep(1, nClasses) # class-specific baseline hazard
  )



  # smart initialization to make separation of class intercepts more likely (intended for 3 classes)
  pars$mu <- c(-1,1) * pars$sigma.b

  currentPars <- unlist(pars)
  
  iter <- 1
  
  mcSamples <- 10
  
  minusTwoLogLikelihood <- NA

  dPredictedList <- list()

  nTimesCriterionMet <- 0
  crit <- Inf

  while (nTimesCriterionMet < 3 && iter <= 100) {
    
    flog.trace(paste0('EM iteration ', iter, ', MC samples: ', mcSamples,', pars = ', paste(format(unlist(pars), digits=0, nsmall=4), collapse = ','),
                      ', crit = ', format(crit, digits = 4)) )

    # Monte Carlo Expectation step
    # Initialize the list of simulated class memberships
    for (mc in 1:mcSamples) {
      dPredictedList[[mc]] <- matrix(NA_integer_, nrow = n, ncol = nClasses)
    }

    #simulate
    for (i in 1:nSubjects) {

      # find conditional distribution of ci for this subject
      pr <- rep(NA_real_, nClasses)
      for (k in 1:nClasses) {
        c1 <- pars$lambda[k] * exp( W[5*(i-1) + 1:5,] %*% pars$alpha + time[5*(i-1) + 1:5] * pars$alpha_time )
        c2 <- - pars$lambda[k] / pars$alpha_time * (exp(pars$alpha_time * time[5*(i-1) + 1:5])-1) * exp(W[5*(i-1) + 1:5,] %*% pars$alpha)

        cik <- rep(0, nClasses)
        cik[k] <- 1

        pr[k] <- dmultinormal(x = y[i, r[i,] == 1 ], #we take only the observed values into account
                              mean = (X[5*(i-1) + 1:5,] %*% pars$beta + as.vector(c(0,pars$mu) %*% cik))[r[i,] == 1],
                              sigma = as.vector( pars$sigma^2 * diag( sum(r[i,]) ) + pars$sigma.b^2)) *
          prod( (c1/(1-pars$theta * c2))^(1-r[i,]) * (1-pars$theta * c2)^(-1/pars$theta) ) *
          dmultinomial( x = cik,
                        size = 1,
                        prob = softmax(c(0,pars$eta)))

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
      
      beta <- x[1:p]
      alpha_time <- x[p + 1]
      alpha <- x[p + 2:pr]
      lsigma.b <- x[p + pr + 1]
      lsigma <- x[p + pr + 2]
      theta <- x[p + pr + 3]
      mu <- x[p + pr + 3 + 1:(nClasses-1)]
      eta <- x[p + pr + 3 + (nClasses-1) + 1:(nClasses-1)]
      llambda <- x[p + pr + 3 + 2*(nClasses-1) + 1:nClasses]

      sigma.b <- exp(lsigma.b)
      sigma <- exp(lsigma)
      lambda <- exp(llambda)

      c1wlk <- exp(alpha_time * time + W %*% alpha)
      c2wlk <- - 1/ alpha_time * (exp(W %*% alpha)-1) * exp(W %*% alpha)

      tmp <- lapply(dPredictedList, 
                    function(dObj) {
                      lambdak <- rep(as.vector(dObj %*% lambda), each = ni)
                      c1 <- c1wlk * lambdak
                      c2 <- c2wlk * lambdak

                      ll <- 0
                      for (i in 1:nSubjects) {
                        if(nio > 0) { #do not execute this code if subject i only has a baseline observation
                          ll <- ll +
                            dmultinormal(x = y[i, r[i,] == 1 ], #we take only the observed values into account
                                         mean = ( X %*% beta + as.vector(c(0,mu) %*% dObj[i,]) )[r[i,] == 1],
                                         sigma = as.vector( sigma^2 * diag( sum(r[i,]) ) + sigma.b^2 ),
                                         log = T)
                        }
                      }
                          ll <- ll + sum( ((1-as.vector(t(r))) * log( c1/(1-theta*c2) ) - (1/theta)*log(1 - theta*c2))) +
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

      beta <- x[1:p]
      alpha_time <- x[p + 1]
      alpha <- x[p + 2:pr]
      lsigma.b <- x[p + pr + 1]
      lsigma <- x[p + pr + 2]
      theta <- x[p + pr + 3]
      mu <- x[p + pr + 3 + 1:(nClasses-1)]
      eta <- x[p + pr + 3 + (nClasses-1) + 1:(nClasses-1)]
      llambda <- x[p + pr + 3 + 2*(nClasses-1) + 1:nClasses]

      sigma.b <- exp(lsigma.b)
      sigma <- exp(lsigma)
      lambda <- exp(llambda)

      out <- lapply(dPredictedList, 
                    function(dObj) {
                      grad <- grad0
                      for (i in 1:nSubjects) {
                        #precalculate quantities
                        nio <- sum(r[i,] == 1)
                        if(nio > 0) { #do not execute this code if subject i only has a baseline observation
                          Xio <- cbind(
                            X[5 * (i - 1) + 1:5,][r[i,] == 1,],
                            matrix(rep(dObj[i, -1], each = nio), nrow = nio)
                          )

                          invSigma <- solve(sigma^2 * diag(nio) + sigma.b^2)
                          ei <- doList[[i]]$y - Xio %*% c(beta, mu)

                          tau <- invSigma %*%
                            (ei %*% t(ei) - (sigma^2 * diag(nio) + sigma.b^2)) %*%
                            invSigma

                          # grad.beta and grad.mu can be calculated at once (score equations as derived by Jennrich and Schluchter)
                          grad[, c(1:p, p + pr + 3 + 1:(nClasses - 1))] <- grad[, c(1:p, p + pr + 3 + 1:(nClasses - 1))] + t(Xio) %*% invSigma %*% ei
                          # then, we have the variance components

                          grad[, p + pr + 1] <- grad[, p + pr + 1] + 1 / 2 *
                            sigma.b *
                            tr(tau %*% (matrix(2 * sigma.b, nrow = nio, ncol = nio)))
                          grad[, p + pr + 2] <- grad[, p + pr + 2] + 1 / 2 * sigma * tr(tau %*% (2 * sigma * diag(nio)))
                        }
                      }

                      #remaining components can be calculated globally:
                      lambdak <- rep(as.vector(dObj %*% lambda), each = ni)

                      #grad.alpha_time
                      grad[p+1] <- -sum(
                        (exp( W %*% alpha ) * (exp( time * alpha_time ) - 1) * ((as.vector(t(r))-1)*theta - 1) * lambdak +
                          time*alpha_time * ((as.vector(t(r))-1)*alpha_time + exp( W %*% alpha ) * (exp(time*alpha_time) + theta - d$r*theta)*lambdak))/
                          (alpha_time*(alpha_time + exp( W %*% alpha )*(exp( time*alpha_time ) - 1)*theta*lambdak)))

                      #grad.alpha_time ... pr
                      for (j in 2:pr) {
                        grad[p+j] <- -sum(
                        ((W[,j]*((as.vector(t(r))-1) * alpha_time + exp(W %*% alpha) * (exp(time*alpha_time) - 1)*lambdak))/
                          (alpha_time + exp( W %*% alpha ) * (exp( time*alpha_time) - 1) * theta * lambdak)
                        ))
                      }

                      #grad.theta
                      grad[p+pr+3] <- sum(
                        ((exp(W %*% alpha) * (exp(time*alpha_time) - 1) * theta * ((as.vector(t(r))-1)*theta -1)*lambdak) /
                          (alpha_time + exp(W %*% alpha) * (exp(time*alpha_time) - 1) * theta * lambdak) +
                          log(1 + (exp(W %*% alpha) * (exp(time*alpha_time) - 1) * theta * lambdak)/alpha_time)
                        )[-seq(from = 1, to = nrow(d), by = ni)])/theta^2

                      #grad.eta
                      grad[,p+pr+3 + (nClasses-1) + 1:(nClasses-1)] <- colSums(dObj[,-1] - rep(softmax(c(0,eta))[-1], each = nrow(dObj)))

                      #grad.lambda
                      dlambdak <- ((1-as.vector(t(r)))*alpha_time - exp(W %*% alpha)*(exp(time*alpha_time)-1)*lambdak ) /
                        (alpha_time + exp(W %*% alpha)*(exp(time*alpha_time)-1)*theta*lambdak )

                      grad[,p+pr+3 + 2*(nClasses-1) + 1:nClasses] <-
                        colSums(matrix(rep(dlambdak, times = nClasses), ncol = nClasses) *
                                  dObj[rep(seq_len(nrow(dObj)), each=ni),])

                      return(grad)
                    }) %>% bind_rows %>% summarize_all(mean) %>% unlist

      -2*out
    }

    ctrl <- list(maxit=100)
    if(nTimesCriterionMet == 2) ctrl$maxit <- 250

    res <- optim(par=unlist(pars, use.names = F), 
                 fn=minusTwoLogLikelihood,
                 gr=minusTwoScore,
                 method = 'BFGS',
                 control = ctrl,
                 hessian = FALSE
    )

    if(res$convergence == 1 && nTimesCriterionMet >= 2) flog.error('class likelihood did not converge within 250 iterations in final step')
    if(res$convergence > 1) flog.error(paste('class likelihood did not converge, code',res$convergence))
    
    pars <- list(beta = res$par[1:p],
                 alpha_time = res$par[p + 1],
                 alpha = res$par[p + 2:pr],
                 sigma.b = exp(res$par[p+pr+1]),
                 sigma = exp(res$par[p+pr+2]),
                 theta = res$par[p+pr+3],
                 mu = res$par[p+pr+3 + 1:(nClasses-1)],
                 eta = res$par[p+pr+3 + (nClasses-1) + 1:(nClasses-1)],
                 lambda = exp(res$par[p+pr+3 + 2*(nClasses-1) + 1:nClasses])
    )

    previousPars <- currentPars
    currentPars <- unlist(pars)

    mcSamples <- floor(min(mcSamples * 1.3, 750)) #increase the sample size slowly
    iter <- iter + 1

    #stopping criteria calculation
    crit <- coalesce(mean( (previousPars-currentPars)^2 )/(mean(previousPars^2)+1e-3), Inf)
    if(crit<1e-4) {
      nTimesCriterionMet <- nTimesCriterionMet + 1
    } else {
      nTimesCriterionMet <- 0
    }
  }
  
  flog.trace(paste0('EM result for class: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ',')))
  
  hess <- hessian(func = minusTwoLogLikelihood,
                  x = unlist(pars, use.names = F))

  return(list(res = res$value, pars = pars, hess = hess))
}

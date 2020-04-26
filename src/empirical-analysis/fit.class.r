#Fit the class model of Lin, McCulloch, and Rosenheck
fit.class <- function(y, r, X, W, nClasses, init, hessMethod = "Richardson") {

  #check inputs
  stopifnot( nClasses >= 2 )
  stopifnot(sum(colnames(W) == "time")==1)

  ni <- ncol(y)
  n <- nrow(y)
  p <- ncol(X)
  pr <- ncol(W)

  time <- W[, colnames(W) == "time"]
  W <- as.matrix(W[, colnames(W) != "time"])

  # prepare grad0 (initial gradient)
  grad0 <- data.frame(0)
  for (j in 1:p) {
    grad0[[paste0('grad.beta',j)]] <- 0
  }
  grad0[['grad.alpha_time']] <- 0
  for (j in 1:(pr-1)) {
    grad0[[paste0('grad.alpha',j)]] <- 0
  }
  grad0[['grad.sigma.b']] <- 0
  grad0[['grad.sigma']] <- 0
  grad0[['grad.theta']] <- 0
  for (k in 2:nClasses) {
    grad0[[paste0('grad.mu',k)]] <- 0
  }
  for (k in 2:nClasses) {
    grad0[[paste0('grad.eta',k)]] <- 0
  }
  for (k in 1:nClasses) {
    grad0[[paste0('grad.lambda',k)]] <- 0
  }
  grad0 <- grad0[,-1]

  pars <- list(beta = init[1:p],
               alpha_time = 1e-2, # time parameter of the hazard function - the likelihood has a removable discontinuity at alpha_time==0 so don't start exactly at it
               alpha = rep(0, pr-1), # remaining hazard parameters
               sigma.b = init[p+1],
               sigma = init[p+2],
               theta = 1, # variance of the frailty distribution
               mu = rep(0, nClasses-1),
               eta = rep(0, nClasses-1),
               lambda = rep(1, nClasses) # class-specific baseline hazard
  )

  # smart initialization to make separation of class intercepts more likely (+1, -1, +2, -2 sigma.b, etc, truncated at number of classes)
  pars$mu <- seq_len(nClasses)[-1] %/% 2 * (-1)^seq_len(nClasses)[-1] * pars$sigma.b

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
    # Initialize the list of simulated datasets
    for (mc in 1:mcSamples) {
      dPredictedList[[mc]] <- matrix(NA_integer_, nrow = n, ncol = nClasses)
    }

    #simulate
    for (i in 1:n) {

      # find conditional distribution of ci for this subject
      prob <- rep(NA_real_, nClasses)
      for (k in 1:nClasses) {
        c1 <- pars$lambda[k] * exp(pars$alpha_time * time[1:ni + ni*(i-1)] + W[1:ni + ni*(i-1),] %*% pars$alpha)
        c2 <- - pars$lambda[k] / pars$alpha_time * (exp(pars$alpha_time * time[1:ni + ni*(i-1)])-1) * exp(W[1:ni + ni*(i-1),] %*% pars$alpha)

        cik <- rep(0, nClasses)
        cik[k] <- 1

        if(sum(r[i,]==1) > 0) {
          prob[k] <- dmultinormal(x = y[i, r[i,] == 1],
                                  mean = (pars$beta[1] +
                                    X[1:ni + ni * (i - 1),] %*% pars$beta +
                                    as.vector(c(0, pars$mu) %*% cik))[r[i,] == 1],
                                  sigma = as.vector(pars$sigma^2 * diag(sum(r[i,] == 1)) + pars$sigma.b^2)) *
            prod(((c1 / (1 - pars$theta * c2))^(1 - r[i,]) * (1 - pars$theta * c2)^(-1 / pars$theta))) *
            dmultinomial(x = cik,
                         size = 1,
                         prob = softmax(c(0, pars$eta)))
        } else {
          # if there are no observed values, there's no normal density
          prob[k] <- prod(((c1 / (1 - pars$theta * c2))^(1 - r[i,]) * (1 - pars$theta * c2)^(-1 / pars$theta))) *
            dmultinomial(x = cik,
                         size = 1,
                         prob = softmax(c(0, pars$eta)))
        }
      }

      ci <- rmultinomial(n = mcSamples, size = 1, prob = prob)

      # Save values
      for (mc in 1:mcSamples) {
        dPredictedList[[mc]][i,] <- ci[mc,]
      }
    }

    flog.trace(paste0('EM iteration ', iter, ' E step completed'))

    # Minimization step
    minusTwoLogLikelihood <- function(x) {
      
      beta <- x[1:p]
      alpha_time <- x[p+1]
      alpha <- x[p+2:pr]
      sigma.b <- x[p+pr+1]
      sigma <- x[p+pr+2]
      theta <- x[p+pr+3]
      mu <- x[p+pr+3 + 1:(nClasses-1)]
      eta <- x[p+pr+3 + (nClasses-1) + 1:(nClasses-1)]
      lambda <- x[p+pr+3 + 2*(nClasses-1) + 1:nClasses]

      c1wlk <- exp(alpha_time * time + W %*% alpha)
      c2wlk <- - 1/ alpha_time * (exp(alpha_time * time)-1) * exp(W %*% alpha)

      tmp <- lapply(dPredictedList, 
                    function(dObj) {
                      lambdak <- rep(as.vector(dObj %*% lambda), each = ni)
                      c1 <- c1wlk * lambdak
                      c2 <- c2wlk * lambdak

                      ll <- 0
                      for (i in 1:n) {
                        if (sum(r[i,]==1) > 0) {
                          ll <- ll +
                            dmultinormal(x = y[i, r[i,] == 1],
                                         mean = (beta[1] +
                                           X[1:ni + ni * (i - 1),] %*% beta +
                                           as.vector(c(0, mu) %*% dObj[i,]))[r[i,] == 1],
                                         sigma = as.vector(sigma^2 * diag(sum(r[i,] == 1)) + sigma.b^2),
                                         log = T)
                        }
                      }
                          ll <- ll + sum( (1-as.vector(t(r))) * log( c1/(1-theta*c2) ) - (1/theta)*log(1 - theta*c2)) +
                          sum(dmultinomial( x = dObj,
                                            size = 1,
                                            prob = softmax(c(0,eta)),
                                            log = T))

                      return(ll)

                    })
      
      out <- -2*mean(unlist(tmp))

      if(any(!is.finite(out))) {
        print(tmp)
        print(out)
        print(list(beta = x[1:p],
                 alpha_time = x[p+1],
                 alpha = x[p+2:pr],
                 sigma.b = x[p+pr+1],
                 sigma = x[p+pr+2],
                 theta = x[p+pr+3],
                 mu = x[p+pr+3 + 1:(nClasses-1)],
                 eta = x[p+pr+3 + (nClasses-1) + 1:(nClasses-1)],
                 lambda = x[p+pr+3 + 2*(nClasses-1) + 1:nClasses]))
        print(minusTwoScore(x))
      }

      out
      
    }

    #TODO: fix analytic gradient. Currently the outputs are incorrect, when compared to a numerical approximation.
    minusTwoScore <- function(x) {

      beta <- x[1:p]
      alpha_time <- x[p+1]
      alpha <- x[p+2:pr]
      sigma.b <- x[p+pr+1]
      sigma <- x[p+pr+2]
      theta <- x[p+pr+3]
      mu <- x[p+pr+3 + 1:(nClasses-1)]
      eta <- x[p+pr+3 + (nClasses-1) + 1:(nClasses-1)]
      lambda <- x[p+pr+3 + 2*(nClasses-1) + 1:nClasses]

      out <- lapply(dPredictedList, 
                    function(dObj) {

                      grad <- grad0
                      for (i in 1:n) {
                        if(sum(r[i,]==1) > 0) {
                          #precalculate quantities
                          nio <- sum(r[i,] == 1)
                          Xi <- cbind(
                            matrix( X[1:ni + ni * (i - 1),][r[i,] == 1,],
                              nrow = nio),
                            matrix( rep(dObj[i, -1], each = nio),
                              nrow = nio))

                          invSigma <- solve(sigma^2 * diag(nio) + sigma.b^2)
                          ei <- y[i, r[i,] == 1] - Xi %*% c(beta, mu)

                          tau <- invSigma %*% (ei %*% t(ei) - (sigma^2 * diag(nio) + sigma.b^2)) %*% invSigma

                          # grad.beta and grad.mu can be calculated at once (score equations as derived by Jennrich and Schluchter)
                          grad[, c(1:p, p + pr + 3 + 1:(nClasses - 1))] <- grad[, c(1:p, p + pr + 3 + 1:(nClasses - 1))] + t(Xi) %*% invSigma %*% ei

                          # then, we have the variance components
                          grad[["grad.sigma.b"]] <- grad[["grad.sigma.b"]] + 1/2 * tr(tau %*% (matrix(2 * sigma.b, nrow = nio, ncol = nio)))
                          grad[["grad.sigma"]] <- grad[["grad.sigma"]] + 1/2 * tr(tau %*% (2 * sigma * diag(nio)))
                        }
                      }

                      #remaining components can be calculated globally:
                      lambdak <- rep(as.vector(dObj %*% lambda), each = ni)

                      #grad.theta
                      grad[["grad.theta"]] <- sum(
                        (exp(W %*% alpha) * (exp(time*alpha_time) - 1) * theta * ((as.vector(t(r))-1)*theta -1)*lambdak) /
                          (alpha_time + exp(W %*% alpha) * (exp(time*alpha_time) - 1) * theta * lambdak) +
                          log(1 + (exp(W %*% alpha) * (exp(time*alpha_time) - 1) * theta * lambdak)/alpha_time)
                        )/theta^2

                      #grad.alpha_time
                      grad[, p+1] <- -sum(
                        (exp(W %*% alpha) * (exp(time*alpha_time) - 1) * ((as.vector(t(r))-1)*theta - 1) * lambdak +
                          time*alpha_time * ((as.vector(t(r))-1)*alpha_time + exp(W %*% alpha) * (exp(time*alpha_time) + theta - as.vector(t(r))*theta)*lambdak))/
                          (alpha_time*(alpha_time + exp(W %*% alpha)*(exp(time*alpha_time) - 1)*theta*lambdak))
                        )

                      #grad.alpha
                      for (j in 1:(pr-1)) {
                        grad[,p+1+j] <- -sum(
                          (W[,j] * ((as.vector(t(r)) - 1) * alpha_time + exp(W %*% alpha) *
                            (exp(time * alpha_time) - 1) * lambdak)) /
                            (alpha_time + exp(W %*% alpha) * (exp(time * alpha_time) - 1) * theta * lambdak))
                      }
                      #grad.eta
                      grad[,p+pr+3 + (nClasses-1) + 1:(nClasses-1)] <- colSums(dObj[,-1] - rep(softmax(c(0,eta))[-1], each = nrow(dObj)))

                      #grad.lambda
                      dlambdak <- ((1-as.vector(t(r)))*alpha_time - exp(W %*% alpha)*(exp(time *alpha_time)-1)*lambdak ) /
                        (lambdak * (alpha_time + exp(W %*% alpha)*(exp(time *alpha_time)-1)*theta*lambdak ))
                      grad[,p+pr+3 + 2*(nClasses-1) + 1:nClasses] <-
                        colSums(matrix(rep(dlambdak, times = nClasses), ncol = nClasses) *
                                  dObj[rep(seq_len(nrow(dObj)), each=ni),])

                      return(grad)
                    }) %>% bind_rows %>% summarize_all(mean) %>% unlist

      print(-2*out)
      g <- grad(minusTwoLogLikelihood, x)
      names(g) <- names(out)
      print(g)

      -2*out
    }
    
    lowerBounds <- rep(-Inf, length(unlist(pars)))
    lowerBounds[c(p+pr+1:3, p+pr+3 + 2*(nClasses-1) + 1:nClasses, length(unlist(pars)))] <- 1e-4 #the three variances, and the lambdas, have to be positive

    ctrl <- list(maxit=min(10*iter,100), trace = 3, REPORT = 1)
    # ctrl <- list(maxit=min(10*iter,100))
    if(nTimesCriterionMet == 2) ctrl$maxit <- 250

    res <- optim(par=unlist(pars, use.names = F), 
                 fn=minusTwoLogLikelihood,
                 # gr=minusTwoScore,
                 method = 'L-BFGS-B',
                 control = ctrl,
                 lower = lowerBounds, 
                 hessian = FALSE
    )

    if(res$convergence == 1 && nTimesCriterionMet >= 2) flog.error('class likelihood did not converge within 250 iterations in final step')
    if(res$convergence > 1) flog.error(paste('class likelihood did not converge, code',res$convergence))
    
    pars <- list(beta = res$par[1:p],
                 alpha_time = res$par[p+1],
                 alpha = res$par[p+2:pr],
                 sigma.b = res$par[p+pr+1],
                 sigma = res$par[p+pr+2],
                 theta = res$par[p+pr+3],
                 mu = res$par[p+pr+3 + 1:(nClasses-1)],
                 eta = res$par[p+pr+3 + (nClasses-1) + 1:(nClasses-1)],
                 lambda = res$par[p+pr+3 + 2*(nClasses-1) + 1:nClasses]
    )

    previousPars <- currentPars
    currentPars <- unlist(pars)

    mcSamples <- floor(min(mcSamples * 1.2 + 2, 500)) #increase the sample size slowly
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

  hess <- hessian(minusTwoLogLikelihood, unlist(pars)) #TODO: switch to jacobian of gradient, once gradient is fixed

  return(list(res = minusTwoLogLikelihood(unlist(pars)), pars = pars, hess = hess))
}
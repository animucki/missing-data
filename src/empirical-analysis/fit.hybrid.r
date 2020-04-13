#Fit the class model of Lin, McCulloch, and Rosenheck, but with Bernoulli missingness
fit.hybrid <- function(y, r, X, W, nClasses, init, hessMethod = "Richardson") {

  ni <- ncol(y)
  n <- nrow(y)
  p <- ncol(X)
  pr <- ncol(W)

  pars <- list(beta = init[1:p],
               alpha = rep(0,pr),
               gamma = 0,
               delta = 0,
               sigma.b = init[p+1],
               sigma = init[p+2],
               mu = rep(0, nClasses-1),
               eta = rep(0, nClasses-1)
  )

  # smart initialization to make separation of class intercepts more likely (intended for 3 classes)
  pars$mu <- c(-1,1) * pars$sigma.b
  
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
                      ', crit = ', format(crit, digits = 4)) )
    
    # Monte Carlo Expectation step
    # Initialize the list of simulated intercepts and class memberships
    for (mc in 1:mcSamples) {
      dPredictedList[[mc]] <- list(cDrawMat=matrix(NA_integer_,
                                                   nrow = n,
                                                   ncol = nClasses),
                                   bDraw=rep(NA_real_, n))
    }

    #simulate
    for (i in 1:n) {

      #simulate bi first
      bVec <- arms(n_samples = mcSamples, lower = -10, upper = 10, metropolis = TRUE,
                   log_pdf = function(bi) {
                     out <- 0
                     for (k in 1:nClasses){
                       cik <- rep.int(0,nClasses)
                       cik[k] <- 1

                       out <- out +
                         prod(dnorm( #outcome
                           x = y[i,],
                           mean = X[5*(i-1) + 1:5,] %*% pars$beta +
                             bi + c(0, pars$mu)[k],
                           sd = pars$sigma), na.rm = T) *
                           prod(dbinom( #missingness
                             x = r[i,],
                             size = 1,
                             prob = plogis( W[5*(i-1) + 1:5,] %*% pars$alpha + pars$gamma * bi + pars$delta * c(0, pars$mu)[k])
                           )) *
                           dnorm( #(normal) random intercept
                             x = bi,
                             sd = pars$sigma.b) *
                           dmultinomial( #class-specific intercept
                             x = cik,
                             size = 1,
                             prob = softmax(c(0, pars$eta)))
                     }

                     return(log(out))
                   })

      for (mc in 1:mcSamples) {
        # Find the conditional class membership distribution
        ciConditional <- rep(NA_real_, nClasses)
        for (k in 1:nClasses) {
          cik <- rep.int(0,nClasses)
          cik[k] <- 1

          ciConditional[k] <-
            prod(dnorm( #outcome
              x = y[i,],
              mean = X[5*(i-1) + 1:5,] %*% pars$beta + bVec[mc] + c(0, pars$mu)[k],
              sd = pars$sigma), na.rm = T) *
              prod(dbinom( #missingness
                x = r[i,],
                size = 1,
                prob = plogis(W[5*(i-1) + 1:5,] %*% pars$alpha + pars$gamma * bVec[mc] + pars$delta * c(0, pars$mu)[k])
              )) *
              dnorm( #(normal) random intercept
                x = bVec[mc],
                sd = pars$sigma.b) *
              dmultinomial( #class-specific intercept
                x = cik,
                size = 1,
                prob = softmax(c(0, pars$eta)))

        }

        # Simulate class and save values
        dPredictedList[[mc]]$cDrawMat[i,] <- rmultinomial(n=1, size=1, prob=ciConditional)
        dPredictedList[[mc]]$bDraw[i] <- bVec[mc]
      }
    }

    flog.trace(paste0('EM iteration ', iter, ' E step completed'))
    
    # Minimization step
    minusTwoLogLikelihood <- function(x) {

      beta <- x[1:p]
      alpha <- x[p + 1:pr]
      gamma <- x[p+pr+1]
      delta <- x[p+pr+2]
      lsigma.b <- x[p+pr+3]
      lsigma <- x[p+pr+4]
      mu <- x[p+pr+4 + 1:(nClasses-1)]
      eta <- x[p+pr+4 + (nClasses-1) + 1:(nClasses-1)]

      sigma.b <- exp(lsigma.b)
      sigma <- exp(lsigma)

      tmp <- lapply(dPredictedList, 
                    function(dObj) {
                      out <- sum(dnorm( #outcome
                        x = as.vector(t(y)),
                        mean = X %*% beta + rep(dObj$bDraw + dObj$cDrawMat %*% c(0, mu), each = ni),
                        sd = sigma,
                        log = T), na.rm = TRUE) +
                        sum(dbinom( #missingness
                          x= as.vector(t(r)),
                          size = 1,
                          prob = plogis( W %*% alpha + rep( gamma * dObj$bDraw + delta * dObj$cDrawMat %*% c(0, mu), each = ni)),
                          log = T))+
                        #second, add up subject-level contributions
                        sum(dnorm( #(normal) random intercept
                          x = dObj$bDraw,
                          sd = sigma.b,
                          log = T)) +
                        sum(dmultinomial( #class-specific intercept
                          x = dObj$cDrawMat,
                          size = 1,
                          prob = softmax(c(0,eta)),
                          log = T))
                      out
                    })
      
      out <- -2*mean(unlist(tmp))

      out
      
    }
    
    minusTwoScore <- function(x) {

      beta <- x[1:p]
      alpha <- x[p + 1:pr]
      gamma <- x[p+pr+1]
      delta <- x[p+pr+2]
      lsigma.b <- x[p+pr+3]
      lsigma <- x[p+pr+4]
      mu <- x[p+pr+4 + 1:(nClasses-1)]
      eta <- x[p+pr+4 + (nClasses-1) + 1:(nClasses-1)]

      sigma.b <- exp(lsigma.b)
      sigma <- exp(lsigma)

      out <- lapply(dPredictedList, 
                    function(dObj) {
                      mcVec <- rep(dObj$cDrawMat %*% c(0, mu), each = ni)
                      normalResidual <- as.vector(t(y)) - X %*% beta - rep(dObj$bDraw, each = ni) - mcVec
                      bernoulliResidual <- as.vector(t(r)) -
                        plogis(W %*% alpha + gamma * rep(dObj$bDraw, each = ni) + delta * mcVec)

                      grad <- data.frame(0)
                      for (j in 1:p) {
                        grad[[paste0('grad.beta',j)]] <- sum( X[,j] * normalResidual, na.rm = T )/sigma^2
                      }
                      for (j in 1:pr) {
                        grad[[paste0('grad.alpha',j)]] <- sum( W[,j] * bernoulliResidual )
                      }

                      grad$grad.gamma <- sum( rep(dObj$bDraw, each = ni) * bernoulliResidual )
                      grad$grad.delta <- sum( mcVec * bernoulliResidual )
                      grad$grad.sigma.b <- sum( (dObj$bDraw^2 - sigma.b^2 )/sigma.b^2)
                      grad$grad.sigma <- sum( normalResidual^2 - sigma^2, na.rm = T )/sigma^2
                      
                      #two separate loops ensure the correct order in the output!
                      for(k in 2:nClasses) {
                        grad[[paste0('grad.mu',k)]] <- sum( normalResidual * rep(dObj$cDrawMat[,k], each = ni), na.rm = T )/sigma^2 +
                          sum( delta * rep(dObj$cDrawMat[,k], each = ni) * bernoulliResidual )
                      }
                      
                      for(k in 2:nClasses) {
                        grad[[paste0('grad.eta',k)]] <- sum( dObj$cDrawMat[,k] - softmax(c(0,eta))[k] )
                      }

                      return(grad[,-1])
                    }) %>% bind_rows %>% summarize_all(mean) %>% unlist

      # print(-2*out/grad(minusTwoLogLikelihood,x))

      -2*out

      
    }

    res <- optim(par=unlist(pars, use.names = F), 
                 fn=minusTwoLogLikelihood,
                 gr=minusTwoScore,
                 method = 'BFGS',
                 hessian = FALSE
    )

    pars <- list(beta = res$par[1:p],
                 alpha = res$par[p + 1:pr],
                 gamma = res$par[p+pr+1],
                 delta = res$par[p+pr+2],
                 sigma.b = exp(res$par[p+pr+3]),
                 sigma = exp(res$par[p+pr+4]),
                 mu = res$par[p+pr+4 + 1:(nClasses-1)],
                 eta = res$par[p+pr+1 + (nClasses-1) + 1:(nClasses-1)])

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
  
  flog.trace(paste0('EM result for spm+class: pars = ', paste(format(unlist(pars), digits=4, nsmall=4), collapse = ',') ) )

  #Calculate the Hessian by calculating the Richardson-extrapolated Jacobian of the gradient
  hess <- jacobian(func = minusTwoScore, x = unlist(pars, use.names = F), method = hessMethod)

  return(list(res = res$value, pars = pars, hess = hess))
}

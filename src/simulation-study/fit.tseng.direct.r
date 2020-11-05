fit.tseng <- function(d) {
  key <- as.integer(d$sample[1])
  set.seed(88810L + key)
  flog.debug(paste0('Fitting Tseng model (direct) to sample ', key))

  #Fit ignorable model to find initial values for parameters

  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  pars <- list(beta = fixef(m),
               alpha = c(0,0,0),
               sigma.b1 = as.data.frame(VarCorr(m))$sdcor[1],
               sigma.b2 = as.data.frame(VarCorr(m))$sdcor[1],
               rho.b = 0,
               sigma = sigma(m)
  )

  zGrid <- mvQuad::createNIGrid(dim = 2, type = 'GHe', level = 3)

  ds <- d %>% dplyr::group_split(subject)

  minusTwoLogLikelihood <- function (x) {
    if(length(x) != 10) flog.error('Incorrect input in log likelihood for parametric model')

    beta <- x[1:3]
    alpha <- x[4:6]
    sigma.b1 <- x[7]
    sigma.b2 <- x[8]
    rho.b <- x[9]
    sigma <- x[10]

    Sigma.b <- matrix(c(sigma.b1^2,
                        sigma.b1 * sigma.b2 * rho.b,
                        sigma.b1 * sigma.b2 * rho.b,
                        sigma.b2^2), ncol = 2)

    llik <- 0

    # make 2d grid
    bGrid <- zGrid
    mvQuad::rescale(bGrid, m=c(0,0), C=Sigma.b)

    for (i in seq_along(ds)) {
      # integrate this subject's likelihood, and accumulate
      Xi <- as.matrix(ds[[i]] %>% dplyr::select(time, treatment))
      llik <- llik + log(mvQuad::quadrature(
        f = function (bi){
          fyiobs <- fri <- fbi <- rep(NA_real_, nrow(bi))

          for (node in seq_len(nrow(bi))) {
            fyiobs[node] <- prod(dnorm(x = ds[[i]]$y,
                                       mean = beta[1] + Xi %*% beta[-1] + bi[node,1]), na.rm = T)
            fri[node] <- prod(dbinom(x = ds[[i]]$r,
                                     size = 1,
                                     prob = plogis(alpha[1] + Xi %*% alpha[-1] + bi[node,2])
            )[-1])
            fbi[node] <- mc2d::dmultinormal(x = bi[node,], mean = c(0,0), sigma = as.vector(Sigma.b))
          }
          return(fyiobs * fri * fbi)
        }, bGrid))
    }

    return(-2*llik)
  }

  res <- optim(par = unlist(pars),
               fn = minusTwoLogLikelihood,
               method = 'L-BFGS-B',
               lower = c(rep(-Inf, 6), 1e-4, 1e-4, -1+1e-4, 1e-4),
               upper = c(rep( Inf, 6),  Inf,  Inf,  1-1e-4, Inf),
               control = list(ndeps = c(rep(1e-3, 6), rep(1e-6, 4)))
  )

  hh <- numDeriv::hessian(func = minusTwoLogLikelihood,
                          x = res$par,
                          method = 'Richardson',
                          method.args = list(eps = 1e-6, d = 1e-6))

  if (all(eigen(hh)$values > 0) & (kappa(hh) < 1e16)) {
    out <- c(key, res$par[c(1,2,3,7,10)], sqrt(diag(2*solve(hh))[c(1,2,3,7,10)]) )
  } else {
    out <- c(key, res$par[c(1,2,3,7,10)], rep(NA_real_, 5) )
  }

  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')

  as.data.frame(as.list(out))
}

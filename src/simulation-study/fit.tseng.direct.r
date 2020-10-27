#Sigma.b <- matrix(c(2,0.9,0.9,2),2)
#bGrid <- mvQuad::createNIGrid(dim = 2, type = 'GHe', level = 10)
#mvQuad::rescale(bGrid, m=c(0,0), C=Sigma.b)
#mvQuad::quadrature(
#  f = function(b) b[,1]^2 * b[,2]^2 * mc2d::dmultinormal(x = b, mean = c(0,0), sigma = as.vector(Sigma.b)),
#  grid = bGrid)
#
#xGrid <- mvQuad::createNIGrid(dim = 1, type = 'GHe', level = 20)
#mvQuad::quadrature(function (x) x^2 * mc2d::dmultinormal(x), xGrid)

fit.tseng <- function(d) {
  key <- as.integer(d$sample[1])
  set.seed(88810L + key)
  flog.debug(paste0('Fitting Tseng model (direct) to sample ', key))

  nTimePoints <- d %>% group_by(subject) %>% summarize(n=n()) %>% pull(n) %>% max

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

    # make grid
    bGrid <- zGrid
    mvQuad::rescale(bGrid, m=c(0,0), C=Sigma.b)

    llik <- 0
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
               control = list(trace=3, REPORT=1)
  )

  hh <- numDeriv::hessian(func = minusTwoLogLikelihood,
                            x = x0)

  out <- c(key, res$par[c(1,2,3,7,10)], sqrt(diag(2*solve(hh))[c(1,2,3,7,10)]) )

  as.data.frame(as.list(out))
}

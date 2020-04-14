#' Fit the Tsonaka et al. semi-parametric model
fit.npsp <- function(y, r, X, W, init) {

  ni <- ncol(y)
  n <- nrow(y)
  p <- ncol(X)
  pr <- ncol(W)

  data. <- y
  miss. <- r

  source("src/tsonaka_etal/SupportFunsIntScaled.R")
  source("src/tsonaka_etal/SupportFunsStdErrorsScaled.R")

  model.1 <- NPSP(data.=y, miss.=r, X=X[,-1], W=W,
                  betas = init[1:p], sigma2 = init[p+1]^2, sigmab = init[p+2]^2,
                  alphas = rep(0,pr), delta = 0,
                  gp = 161, bound = 3, reltol = 1e-08, epsilon = 0.001, iter = 1000)

  # Standard error calculation
  NPbetas <- model.1$bts
  NPalphas <- model.1$alps
  NPalphas2 <- model.1$alphas
  NPgamma <- model.1$gamma.
  NPvb <- model.1$vb
  NPmb <- model.1$mb
  NPsig <- sqrt(model.1$sig)
  NPmiss <- model.1$per.mis
  NPsb <- model.1$sigmab

  Sup <- model.1$sup.points
  pi. <- Sup[, 2]
  gridp <- Sup[, 1] - NPmb

  # Profile
  y <- c(t(data.))
  n <- nrow(data.) # sample size
  p <- ncol(data.)
  na <- ncol(W)
  gp <- length(gridp)

  #very hacky fix to "object X not found" error
  X <- X[,-1]
  environment(min.logALLProfile) <- environment(grad.ALLProfile) <- environment(hess.parsProfile) <- environment()

  hess <- hess.parsProfile(c(c(NPbetas, NPsig^2, NPalphas, NPgamma)))

  return(list(res = -2*model.1$logLik., pars = c(beta = NPbetas, sigma.b = sqrt(NPsb), sigma = NPsig, alpha = NPalphas, gamma = NPgamma), hess = hess))
}

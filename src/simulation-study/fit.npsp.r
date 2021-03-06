#' Fit the Tsonaka et al. semi-parametric model
fit.npsp <- function(d) {

  key <- as.integer(d$sample[1])
  flog.debug(paste0('Fitting semi-parametric random effect model to sample ', key))

  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)

  ni <- d %>% group_by(subject) %>% summarize(n=n()) %>% pull(n) %>% max

  data. <- matrix(d$y, byrow = TRUE, ncol = ni)
  miss. <- matrix(d$r, byrow = TRUE, ncol = ni)[,-1]

  X <- as.matrix(d[,c("time","treatment")])
  W <- cbind(1, as.matrix(d[,c("time","treatment")]))[-seq(from=1,to=nrow(d),by=ni),]

  source("src/tsonaka_etal/SupportFunsIntScaled.R")
  source("src/tsonaka_etal/SupportFunsStdErrorsScaled.R")

  model.1 <- NPSP(data.=data., miss.=miss., X=X, W=W,
                  betas = fixef(m), sigma2 = sigma(m)^2, sigmab = as.data.frame(VarCorr(m))$vcov[1],
                  alphas = c(0,0,0), delta = 0,
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
  environment(min.logALLProfile) <- environment(grad.ALLProfile) <- environment(hess.parsProfile) <- environment()

  hess <- hess.parsProfile(c(c(NPbetas, NPsig^2, NPalphas, NPgamma)))
  se <- sqrt(diag(2*solve(hess))[1:4])

  out <- c(key
    , model.1$bts, sqrt(model.1$sigmab), sqrt(model.1$sig), se[1:3], NA, se[4] )
  names(out) <- c('sample','intercept', 'time', 'treatment', 'sigma.b', 'sigma',
                  'se.intercept', 'se.time', 'se.treatment', 'se.sigma.b', 'se.sigma')
  as.data.frame(as.list(out))
}

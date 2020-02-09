fit.ignorable <- function(d) {
  key <- d$sample[1]
  flog.debug(paste0('Fitting ignorable model to sample ', key))
  
  m <- lmer(y ~ (1|subject) + time + treatment,
            data=d, REML=F)
  beta <- fixef(m)
  names(beta)[1] <- 'intercept'
  
  seFixed <- diag(vcov(m))
  names(seFixed) <- paste0('se.',names(beta))
  
  #extract Wald standard errors for variance components (not included in output by default)
  dd.ML <- lme4:::devfun2(m, useSc=TRUE, signames=FALSE)
  vv <- as.data.frame(VarCorr(m))
  pars <- vv$sdcor
  names(pars) <- c('sigma.b','sigma')
  hh <- hessian(dd.ML, pars)
  vv <- sqrt(diag(2*solve(hh)))
  names(vv) <- c('se.sigma.b','se.sigma')
  
  as.data.frame(as.list(c(beta, pars, seFixed, vv)))
}
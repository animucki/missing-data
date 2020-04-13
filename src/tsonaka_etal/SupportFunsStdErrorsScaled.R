##### Functions to compute profile likelihood hessian evaluated at beta0 + E(b), grid = grid - E(b), a0 = a0 + E(b)

# minus - log lik as a function of betas, sigma2, alphas, gamma, pi.
min.logALLProfile <- function(pars){
    betas.0 <- pars[1]
    betas <- pars[1:2 + 1]
    sigma2 <- pars[3 + 1]
    alphas <- pars[4:6 + 1]
    delta <- pars[7 + 1]
    fixY <- outer(betas.0 + drop(X %*% betas), gridp, "+") 
    prY <- dnorm(y, mean = fixY, sd = sqrt(sigma2), log = TRUE) 
    indx <- rep(1:n, each = p) 
    Y.b <- rowsum(prY, group = indx, na.rm = TRUE)
    
    fixR <- drop(W %*% alphas)
    R.b <- matrix(0, nrow = n, ncol = gp)
    for(i in 1:gp){
        predR <- matrix(fixR + delta * gridp[i], byrow = TRUE, ncol = p - 1)
        R.b[, i] <- rowSums(miss. * predR + log(1 - plogis(predR)), na.rm = TRUE)
        }
    min.log <- -sum(log(exp(Y.b + R.b) %*% pi.))
    min.log
}

#lg <- min.logALL(c(model.1$bts[-1], model.1$sig, model.1$alphas, model.1$gamma., model.1$sup.points[-1, 2]))
#lg


# minus gradient of the logLik wrt both alphas and delta - CD
grad.ALLProfile <- function(params){            
            cd(params, f = min.logALLProfile)
}
#grALL <- grad.ALL(c(model.1$bts[-1], model.1$sig, model.1$alphas, model.1$gamma., model.1$sup.points[-1, 2]))
#grALL


# hessian wrt both betas and sigma2 - CD
hess.parsProfile <- function(params){            
            cd.vec(params, f = grad.ALLProfile)
}


##### Functions to compute hessian using re-distributing weights evaluated at beta0 + E(b), grid = grid - E(b), a0 = a0 + E(b)

####
# minus - log lik as a function of betas, sigma2, alphas, gamma, pi.
min.logALLThres <- function(pars){
    betas.0 <- pars[1]
    betas <- pars[1:2 + 1]
    sigma2 <- pars[3 + 1]
    alphas <- pars[4:6 + 1]
    delta <- pars[7 + 1]
    pi. <- pars[-c(1:(7 + 1))]
    pi. <- c((1 - sum(pi.)), pi.)
    fixY <- outer(betas.0 + drop(X %*% betas), gridp, "+") 
    prY <- dnorm(y, mean = fixY, sd = sqrt(sigma2), log = TRUE) 
    indx <- rep(1:n, each = p) 
    Y.b <- rowsum(prY, group = indx, na.rm = TRUE)
    
    fixR <- drop(W %*% alphas)
    R.b <- matrix(0, nrow = n, ncol = gp)
    for(i in 1:gp){
        predR <- matrix(fixR + delta * gridp[i], byrow = TRUE, ncol = p - 1)
        R.b[, i] <- rowSums(miss. * predR + log(1 - plogis(predR)), na.rm = TRUE)
        }
    min.log <- -sum(log(exp(Y.b + R.b) %*% pi.))
    min.log
}

#lg <- min.logALL(c(model.1$bts[-1], model.1$sig, model.1$alphas, model.1$gamma., model.1$sup.points[-1, 2]))
#lg


# minus gradient of the logLik wrt both alphas and delta - CD
grad.ALLThres <- function(params){            
            cd(params, f = min.logALLThres)
}
#grALL <- grad.ALL(c(model.1$bts[-1], model.1$sig, model.1$alphas, model.1$gamma., model.1$sup.points[-1, 2]))
#grALL


# hessian wrt both betas and sigma2 - CD
hess.parsThres <- function(params){            
            cd.vec(params, f = grad.ALLThres)
}



### Functions to compute the Hessian for the common SPM

min.logALL <- function(pars){
    betas <- pars[1:3]
    sigma <- pars[4]
    alphas <- pars[5:7]
    gama <- pars[8]
    sigma.b <- pars[9]
    
    # create required data matrices
    miss. <- t(apply(dat, 1, function(x){ x <- ifelse(!is.na(x), 1, 0); x})) # missingness matrix
    data. <- dat
    p <- ncol(data.)
    gfr. <- data.frame(GFR = c(t(data.)), treat = d.[, 3], time. = d.[, 4], patid = d.[, 5])
    gfr <- na.exclude(gfr.)
    Rn <- c(t(miss.)); Rn <- Rn[-seq(from = 1, to = length(Rn), by = p)]
    Rn. <- 1 - Rn
    id <- gfr$patid
    y <- gfr$GFR
    X <- model.matrix( ~ time. * treat - treat, data = gfr)
    attr(X, "assign") <- NULL

    namsX <- colnames(X)
    Z <- X[, c("(Intercept)"), drop = FALSE]
    n <- length(unique(id))
    id <- rep(1:n, rle(id)$length)
    N <- nrow(X)

    W <- model.matrix(~ treat + time., data = gfr.); attr(W, "assign") <- NULL
    W <- W[-seq(from = 1, to = p * nrow(dat), by = p), ]
    namsW <- colnames(W)
    q <- 1 # dim rand-eff
    qy <- 1

    tXX <- crossprod(X)
    tZZ <- lapply(split(Z, id), function(x) crossprod(matrix(x, ncol = qy))); names(tZZ) <- NULL
    tZZ <- matrix(unlist(tZZ), n, qy * qy, TRUE)

    k <- qp
    GH <- gauher(k)
    b <- b. <- as.matrix(expand.grid(lapply(1:q, function(k, u) u$x, u = GH)))
    b <- sqrt(2) * b
    k <- nrow(b)
    bqy <- b[, 1:qy]
    bqy2 <- if(qy == 1) bqy * bqy else t(apply(bqy, 1, function(x) x %o% x))
    Z.tbqy <- Z %*% t(bqy)
    mat.bq <- matrix(rep(b[, q], p - 1), n * (p - 1), k, TRUE)
    wGH <- as.matrix(expand.grid(lapply(1:q, function(k, u) u$w, u = GH)))
    wGH <- sqrt(2) * apply(wGH, 1, prod) * exp(rowSums(b. * b.))

    sigma.y <- sigma
    gammas <- alphas
    alpha <- gama

    mu.y <- c(X %*% betas) + Z.tbqy
    log.p.yb <- rowsum(dnorm(y, mu.y, sigma.y, log = TRUE), id)
    mu.t <- plogis(c(W %*% gammas) + alpha * mat.bq)
    mr <- Rn * log(mu.t) + Rn. * log(1 - mu.t)
    prs. <- array(mr, c(p - 1, n, k))
    log.p.tb <- colSums(prs., 2)
    p.ytb <- exp(log.p.yb + log.p.tb)
    p.b <- dnorm(b, sd = sigma.b)
    wGHn <- wGH * p.b
    p.yt <- c(p.ytb %*% wGHn)
    min.log <- -sum(log(p.yt))
    min.log
}


# minus gradient of the logLik wrt both Y, R - CD
grad.ALL <- function(pars){            
            cd(pars, f = min.logALL)
}
#grALL <- grad.ALL(c(model.1$bts[-1], model.1$sig, model.1$alphas, model.1$gamma., model.1$sup.points[-1, 2]))
#grALL


# hessian wrt both betas and sigma2 - CD
hess.pars <- function(pars){            
            cd.vec(pars, f = grad.ALL)
}

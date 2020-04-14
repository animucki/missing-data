########################################################################################
# program : SupportFuns.R                                                               
# project : Shared Parameter models with Flexible Mixing distributions (NP)             
# purpose : Functions required to fit SPSP model using the VEM and the SP model         
#           using the EM algorithm (g = 0 is also considered for both models)           
# creation date : 07MARCH2007                                                           
# update :                                                                              
# Note: Data are simulated with zero mean for b; distribution assumed for b are:        
#       normal, mixture, log-normal                                                     
# programmer : R. Tsonaka                                                               
########################################################################################


#########################################################################################
# Details for the functions for NPSPM                                                    
# Aim:  Functions to run the SPM for continuous normal data using the VEM (b unspecified)
#       To make the RE to have zero mean we fix the intercept betas.0 at a value and     
#       update it at the end. We do not centralize the grid after each iteration because 
#       the likehood fails to increase. The reason for this is that after the            
#       transformation  from b to b* (with mean 0)implies that the RE density must be    
#       estimated at f(b + mean(b)) but in these b + mean(b) support points the weigts   
#       are not known only the weights for b are known (check in Tsiatis book Chapter 1).
#       The step length is set to a = 1 as long as likelihood increases otherwise is     
#       estimated. Random intercepts only.                                               
# Date: 6/2/2007                                                                         




# Plot of the fitted random effects distribution - lines - NPSP
fitted.b.plotL <- function(prec, sb, sup, D.Y, typ = "l"){
            
    #plot(sup[, 1], sup[, 2], col = 1, cex.lab = 2, cex.axis = 1.9, type = typ, xlab = "", ylab = "")
    #mn. <- sum(sup[, 1] * sup[, 1])
    segments(x0 = sup[, 1], y0 = rep(0, length(sup[, 1])), 
             x1 = sup[, 1], y1 = sup[, 2])
    #mtext("Grid Points", side = 1, line = 4, cex = 1.9)
    #mtext("Density", side = 2, line = 4, cex = 1.9)

}

# Plot of the simulated normal random effects distribution SPM - lines
sim.b.plotNl <- function(prec, sb, n, D.Y){
    x <- seq(-sb, sb, length = prec)
    grid. <- as.matrix(x)
    f <- function(x){
        mu. <- 0
        dens <- dnorm(x, mu., D.Y)
    }
    z1 <- apply(grid., 1, function(x) f(x))
    dim(z1) <- c(prec)
    lines(x, z1, col = 1, cex = 0.45)
}



# summing the elements of a list
matSums <- function(lis){
    ress <- array(data = 0.0, dim = dim(lis[[1]]))
    for(i in seq(along = lis)){
        ress <- ress + lis[[i]]
        }
        ress
        }


# Numerical derivative - central difference - for scalars
cd <- function (x, f, ..., eps = 1e-04) {
    n <- length(x)
    res <- numeric(n)
    ex <- pmax(abs(x), 1)
    for (i in 1:n) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- 2 * max(abs(c(x1[i] - x[i], x2[i] - x[i])))
        res[i] <- diff.f / diff.x
    }
    res
}


# numerical approximation for score vectors - for vectors
cd.vec <- function (x, f, ..., eps = 1e-04) {
    n <- length(x)
    res <- matrix(0, n, n)
    ex <- pmax(abs(x), 1)
    for (i in 1:n) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- 2 * max(abs(c(x1[i] - x[i], x2[i] - x[i])))
        res[, i] <- diff.f / diff.x
    }
    res
} 


# minus log-likelihood
min.log <- function(betas, sigma2, alphas, delta){
        
    fixY <- outer(betas.0 + drop(X %*% betas), gridp, "+") 
    prY <- dnorm(y, mean = fixY, sd = sqrt(sigma2), log = TRUE) 
    indx <- rep(1:n, each = p) 
    Y.b <- rowsum(prY, group = indx, na.rm = TRUE)
    
    fixR <- drop(W %*% alphas)
    R.b <- matrix(0, nrow = n, ncol = gp)
    for(i in 1:gp){
        predR <- matrix(fixR + delta * gridp[i], byrow = TRUE, ncol = ncol(miss.))
        R.b[, i] <- rowSums(miss. * predR + log(1 - plogis(predR)), na.rm = TRUE)
        }
    min.log <- -sum(log(exp(Y.b + R.b) %*% pi.))
    list(min.log = min.log, Y.b = exp(Y.b), R.b = exp(R.b))
}
 
#min.log(betas, sigma2, alphas, delta)$min.log


# minus log-likelihood as a function of both betas and sigma2
min.logbm <- function(params, R.b){
    
    betas <- params[1:nb]
    sigma2 <- exp(params[nb + 1])
    fixY <- outer(betas.0 + drop(X %*% betas), gridp, "+")
    prY <- dnorm(y, mean = fixY, sd = sqrt(sigma2), log = TRUE)
    indx <- rep(1:n, each = p) 
    Y.b <- rowsum(prY, group = indx, na.rm = TRUE)
    min.log <- -sum(log((exp(Y.b) * R.b) %*% pi.))
    min.log
} 
#aa <- min.logbm(c(betas, log(sigma2)), R.b)


# minus gradient of the logLik wrt both betas and sigma2 - CD
grad.bm2 <- function(params, R.b){            
            cd(params, f = min.logbm, R.b = R.b)
}
#grad.bm2(c(betas, log(sigma2)), R.b)


# hessian wrt both betas and sigma2 - CD
hess.bm2 <- function(params, R.b){            
            cd.vec(params, f = grad.bm2, R.b = R.b)
}
#hess.bm2(c(betas, log(sigma2)), alphas, delta, nb, na, X, gridp, y, n, p, W, gp, pi.)


# minus gradient of the logLik wrt betas and sigma2 - CD
grad.bm <- function(params, R.b){
            
            betas <- params[1:nb]
            sigma2 <- exp(params[nb + 1])
            fixY <- outer(betas.0 + drop(X %*% betas), gridp, "+") 
            prY <- dnorm(y, mean = fixY, sd = sqrt(sigma2), log = TRUE) 
            indx <- rep(1:n, each = p) 
            Y.b <- exp(rowsum(prY, group = indx, na.rm = TRUE))
            fyr.b <- Y.b * R.b
            fyrbp <- (Y.b * R.b) * matrix(pi., nrow = n, ncol = length(pi.), byrow = TRUE)
            fyr <- fyr.b %*% pi.
            
            dl.b <- numeric(nb)
            f.1 <- outer(y - betas.0 - drop(X %*% betas), gridp, "-")
            for(bb in 1:nb){
                dl <- f.1 * X[, bb]
                indx <- rep(1:n, each = p) 
                dl <- rowsum(dl, group = indx, na.rm = TRUE)
                dl.b[bb] <- (1 / sigma2) * sum((fyrbp / c(fyr)) * dl)
            }
            dl.s <- rowsum(-1 / 2 + (1 / (2 * sigma2)) * f.1^2, group = indx, na.rm = TRUE)
            
            -c(dl.b, sum((fyrbp / c(fyr)) * dl.s))
}
# compare gradients:
# grad.bm(c(betas, log(sigma2)), R.b)
# grad.bm2(c(betas, log(sigma2)), R.b)



# minus log-likelihood as a function of betas only
min.logbm1 <- function(params){

    betas <- params
    fixY <- outer(betas.0 + drop(X %*% betas), gridp, "+")
    prY <- dnorm(y, mean = fixY, sd = sqrt(sigma2), log = TRUE)
    indx <- rep(1:n, each = p) 
    Y.b <- rowsum(prY, group = indx, na.rm = TRUE)
    min.log <- -sum(log((exp(Y.b) * R.b) %*% pi.))
    min.log
} 
#aa <- min.logbm1(betas)


# minus gradient of the logLik wrt betas only - CD
grad.bm1 <- function(params){            
            cd(params, f = min.logbm1)
}
#gr.bm <- grad.bm1(c(betas), sigma2, y, X, gridp, pi., p, fyr, fyrbp, nb, n, R.b)


# function to be maximized to find an update for a
lgL.fun <- function(alpha){
                pi.n <- pi.
                pi.n[ind.min] <- (1 - alpha) * pi.[ind.min]
                pi.n[ind.max] <- alpha * pi.[ind.min] + pi.[ind.max] 
                fyr <- (Y.b * R.b) %*% pi.n
                sum(log(fyr)) - logLik.
                }

# minus logLik wrt alphas, delta
min.logR <- function(pars){
    delta <- pars[1]
    alphas <- pars[-1]
    fixR <- W %*% alphas
    R.b <- matrix(0, nrow = n, ncol = gp)
    for(i in 1:gp){
        predR <- matrix(fixR + delta * gridp[i], byrow = TRUE, ncol = ncol(miss.))
        R.b[, i] <- exp(rowSums(miss. * predR + log(1 - plogis(predR))))
    }
    -sum(log((Y.b * R.b) %*% pi.))
    
} 

#min.logR(c(delta, alphas), W, n, gp, grid.new, p, miss., Y.b, pi., na)

# minus logLik wrt alphas, delta
grad.R <- function(pars){

    delta <- pars[1]
    alphas <- pars[-1]
    fixR <- W %*% alphas
    indx <- rep(1:n, each = ncol(miss.)) 
    R.b <- matrix(0, nrow = n, ncol = gp)
    Rpw <- matrix(0, nrow = n, ncol = gp)
    Rpm <- matrix(0, nrow = n, ncol = gp)
    dl.b <- numeric(na)  
    for(bb in 1:na){
    for(i in 1:gp){
        predR <- matrix(fixR + delta * gridp[i], byrow = TRUE, ncol = ncol(miss.))
        Rpw[, i] <- rowsum(c(t((miss. - plogis(predR)))) * W[, bb], group = indx, na.rm = TRUE)
        Rpm[, i] <- rowsum(c(t((miss. - plogis(predR)))) * gridp[i], group = indx, na.rm = TRUE)
        R.b[, i] <- exp(rowSums(miss. * predR + log(1 - plogis(predR))))
        }
    fyr.b <- Y.b * R.b
    fyrbp <- (Y.b * R.b) * matrix(pi., nrow = n, ncol = length(pi.), byrow = TRUE)
    fyr <- fyr.b %*% pi.
    dl.b[bb] <- sum((fyrbp / c(fyr)) * Rpw)   
    }
    dl.d <- sum((fyrbp / c(fyr)) * Rpm) 
    
    -c(dl.d, dl.b)
} 

#grad.R(c(delta, alphas), W, n, gp, grid.new, p, miss., Y.b, pi., na)


# minus gradient of the logLik wrt alphas, delta
grad.R2 <- function(pars){
                cd(pars, f = min.logR)
                }   
#grad.R2(c(delta, alphas), W, n, gp, grid.new, p, miss., Y.b, pi.)

# second derivative of the logLik wrt alphas, delta
hess.R <- function(pars){
                cd.vec(pars, f = grad.R)
                }   


# Main function
NPSP <- function(data., miss., X, W, betas, sigma2, sigmab, alphas, delta, gp, bound, reltol = 1e-08, epsilon = 0.001, iter = 1000){
    
    nb <- ncol(X)
    y <- c(t(data.))
    n <- nrow(data.) # sample size
    p <- ncol(data.)
    na <- ncol(W)
    gridp <- sigmab * seq(-bound, bound, length = gp) # grid
    pi. <- rep(1/gp, gp) # weights
    # measurement
    betas.0 <- betas[1]
    betas <- betas[-1]
    conv <- FALSE
    grid.old <- gridp
    pi.old <- numeric(length(pi.))  
    #grid.new <- gridp
    env <- environment()
    environment(min.log) <- environment(min.logbm) <- environment(grad.bm) <- env
    environment(min.logR) <- environment(grad.R) <- environment(grad.R2) <- environment(lgL.fun) <- env
    
    for(it in 1:iter){    
        
        likEval <- min.log(betas, sigma2, alphas, delta)
        Y.b <- likEval$Y.b
        R.b <- likEval$R.b
        #used only to return results
        fyr.b <- Y.b * R.b
        fyrbp <- (Y.b * R.b) * matrix(pi., nrow = n, ncol = length(pi.), byrow = TRUE)
        fyr <- fyr.b %*% pi.

        if(it > 1){
            logLik.. <- - likEval$min.log
            difL <- logLik.. - logLik.
            loglikP <- logLik. 
            }
                    
        if(it > 1 && (difL < 0 & abs(difL) > reltol * (abs(logLik.) + reltol))){
            # cat("\niteration:", it, "logLik:", logLik..)
            stop("\n************Likelihood failed to increase**************\n") 
            #warning("\nLikelihood failed to increase!!!!!!!!!!!!!!\n")
            }

        if(!conv){
            fyr <- (Y.b * R.b) %*% pi.         
            gf.val <- colMeans(Y.b * R.b / c(fyr))
            ind.max <- which(gf.val == max(gf.val))[1] #; cat("\nmax", ind.max)
            ind.min <- which(gf.val == min(gf.val))[1] #; cat("\nmin", ind.min)
        }
        logLik. <- - likEval$min.log
        # cat("\niteration:", it, "logLik:", logLik.)
        # if(it > 1) cat("\nDif:", difL)
        # cat("\ngradient", gf.val[ind.max])
        if(gf.val[ind.max] <= 1 + epsilon){
            # cat("\nConverged!")
            conv <- TRUE
            }
        #if(it > 1 && (conv & difL < 1e-7)){
        if(it > 1 && (conv & abs(difL) < reltol * (abs(logLik.) + reltol))){
        
            # cat("\nConverged!")
            break
            }
        if(!conv){
            a. <- if(lgL.fun(1) > 0) 1 else{
                     optimize(f = lgL.fun, interval = c(0, 1), maximum = TRUE)$maximum}
            # cat("\nIncr:", lgL.fun(1))
            # cat("\nIncr:", lgL.fun(a.))
            # cat("\nalpha:", a.)

            pi.n <- pi.
            pi.n[ind.min] <- (1 - a.) * pi.[ind.min]

            pi.n[ind.max] <- c(a. * pi.[ind.min]) + c(pi.[ind.max]) 
            pi. <- pi.n
            # cat(pi.n[ind.min], pi.n[ind.max], "\n")
            
            ind.zero <- which(pi. < sqrt(.Machine$double.eps))
            # ind.zero <- which(pi. < 1e-06)
            
            if(length(ind.zero !=0 )){
                pi.[ind.zero] <- 0
                pi. <- pi./sum(pi.)
                pi. <- pi.[-ind.zero]
                gp <- gp - length(ind.zero)
                # cat("\nPoints:", gp)
                # print(ind.zero)
                gridp <- gridp[-ind.zero]
            }

            id.grid <- which(grid.old %in% gridp)
            pi.old[id.grid] <- pi.  
            # cat("E(b):", gridp %*% pi., "\n")
            # cat("\nWeights:", pi.)
            # cat("\nSum:", sum(pi.))
            }
        # cat("Grid:", gridp, "\n")
        # grid.new <- gridp + gridp %*% pi.
        # cat("E(b.new):", grid.new %*% pi., "\n")
###
        likEval <- min.log(betas, sigma2, alphas, delta)
        Y.b <- likEval$Y.b
        R.b <- likEval$R.b
        
# if(-likEval$min.log > logLik.) cat("\n VEM increased the likelihood:", -likEval$min.log) else cat("\n ***** VEM DID NOT increase the likelihood ******", -likEval$min.log)
               

# update betas + sigma2
bts <- optim(par = c(betas, log(sigma2)), fn = min.logbm, gr = grad.bm, 
             method = "BFGS", control = list(maxit = 5, 
             parscale = c(rep(1, length(betas)), 0.1)), R.b = R.b)

#bl.1 <- -min.logbm(c(betas, log(sigma2)), R.b = R.b)
#bl.2 <- -min.logbm(bts$par, R.b = R.b)
#cat("\nb start:", bl.1)
#cat("\nb final:", bl.2)

#if(bl.2 > bl.1) cat("\nbs increased the likelihood") else cat("\n***** bs DID NOT increase the likelihood ******")

betas <- bts$par[1:nb]
sigma2 <- exp(bts$par[nb + 1])
# cat("\nbetas:", betas)
# cat("\nsigma2:", sigma2)
# cat("\nb conv:", bts$conv)
# cat("\nb conv:", bts$counts)

likEval <- min.log(betas, sigma2, alphas, delta)
Y.b <- likEval$Y.b
R.b <- likEval$R.b
# cat("\nLog-Lik after bs:", -likEval$min.log)


# update alphas + delta
alps <- optim(par = c(delta, alphas), fn = min.logR, gr = grad.R, method = "BFGS", 
              control = list(maxit = 5, parscale = rep(0.1, length(c(delta, alphas)))))

#al.1 <- -min.logR(c(delta, alphas))
#al.2 <- -min.logR(alps$par)
#cat("\na start:", al.1)
#cat("\na final:", al.2)
#if(al.2 > al.1) cat("\nas increased the likelihood") else cat("\n***** as DID NOT increase the likelihood ******")
# cat("\nalphas conv:", alps$conv)
# cat("\nalphas counts:", alps$counts)

alphas <- alps$par[-1]
delta <- alps$par[1]
# cat("\nalphas:", alphas)
# cat("\ndelta:", delta)

}

list(sup.points = cbind(params = gridp, weight = pi., gr = gf.val), 
     iter = it, 
     sup.points.all = cbind(params = grid.old, weight = pi.old),
     max.gr = max(gf.val[pi. > sqrt(.Machine$double.eps)]), 
     bts = c(betas.0 + sum(gridp * pi.), betas),  
     alps = c(alphas[1] + delta * sum(gridp * pi.), alphas[-1]),
     gamma. = delta, mb = sum(gridp * pi.), 
     vb = sum((gridp - sum(gridp * pi.))^2 * pi.), 
     logLik. = logLik., 
     post = cbind(ID = 1:n, Class = unlist(apply(c(1/fyr) * (fyrbp), 1, function(x) (gridp)[which(x %in% max(x))]))),
     gp = gp, y = y, X = X, W = W, Y.b = Y.b, R.b = R.b,
     sigmab = sigmab,
     alphas = alphas,
     betas.0 = betas.0,
     fyr = fyr, 
     fyrbp = fyrbp,
     sig = sigma2, 
     data. = data., 
     # sim.b = summary(simul$b),
     per.mis = sum(is.na(data.)) / (n * p)
     # bouns = max(simul$b)
)

}

#mod <- NPSP(data. = simul$dat, miss. = miss, X = X, W = W, betas = cfs.lme, sigma2 = sigma^2, sigmab = sigmab, 
#            alphas = cfs.lmer, delta = 0.7, gp = 61, bound = 7, reltol = 1e-08, epsilon = 0.001, iter = 1000)
            
#plot(mod$sup.points[, 1], mod$sup.points[, 2], xlim = c(-5, 10), ylim = c(0, 1), type = "n")
#segments(x0 = mod$sup.points[, 1], y0 = rep(0, length(mod$sup.points[, 1])), 
#         x1 =mod$sup.points[, 1],  y1 = mod$sup.points[, 2], xlim = c(-5, 10), ylim = c(0, 1))

#lines(seq(-1, 5, by=0.1)-1.377128, dlnorm(seq(-1, 5, by=0.1), 0, 0.8))


#exp(2*in.pars$loc + in.pars$scale.^2)*(exp(in.pars$scale.^2) - 1)
#exp(2*0 + 0.97^2)*(exp(0.97^2) - 1)





# Gauss-Hermite: returns x, w so that \int_{-\infty}^{\infty} exp(-x^2) f(x) dx \doteq \sum w_i f(x_i)
gauher <- function(n){ 
            EPS <- 3.e-14
            PIM4 <- .7511255444649425
            MAXIT <- 10
            m <- trunc((n + 1) / 2)
            x <- w <- rep(-1, n)
            for (i in 1:m) {
                if (i == 1) {
                z <- sqrt(2 * n + 1) - 1.85575 * (2 * n + 1)^(-.16667)
                } else if(i == 2) {
                z <- z - 1.14 * n^.426 / z    
                } else if (i == 3) {
                z <- 1.86 * z - .86 * x[1]
                } else if (i == 4) {
                z <- 1.91 * z - .91 * x[2]
                } else {
                z <- 2. * z - x[i - 2]
                }
            for (its in 1:MAXIT){
                    p1 <- PIM4
                    p2 <- 0.
                    for (j in 1:n){
                            p3 <- p2
                            p2 <- p1
                            p1 <- z * sqrt(2. / j) * p2 - sqrt((j - 1) / j) * p3
                            }
                            pp <- sqrt(2. * n) * p2 
                            z1 <- z
                            z <- z1 - p1 / pp
                            if(abs(z - z1) <= EPS) break
                            }
                    x[i] <- z
                    x[n + 1 - i] <- -z
                    w[i] <- 2 / (pp * pp)
                    w[n + 1 - i] <- w[i]
                    }
            list(x = x, w = w)
            }

# Support function for the SP model

opt.random <- function(thetas){
    log.p.b <- dnorm(b, sd = exp(thetas), log = TRUE)
    -sum(c(p.byt %*% (log.p.b * wGHn)))
}

opt.missin <- function(thetas){
    gammas <- thetas
    mu.t <- plogis(c(W %*% gammas) + alpha * mat.bq)
    mr <- Rn * log(mu.t) + Rn. * log(1 - mu.t)
    prs. <- array(mr, c(ncol(miss.), n, k))
    log.p.tb <- colSums(prs., 2)
    p.bytn <- p.byt * log.p.tb
    -sum(c(p.bytn %*% wGHn))
}

opt.missin2 <- function(alph){
    mu.t <- plogis(c(W %*% gammas) + alph * mat.bq)
    mr <- Rn * log(mu.t) + Rn. * log(1 - mu.t)
    prs. <- array(mr, c(ncol(miss.), n, k))
    log.p.tb <- colSums(prs., 2)
    p.bytn <- p.byt * log.p.tb
    -sum(c(p.bytn %*% wGHn))
}
            
# Function to fit the SP model with normal RE
# sigma is sd
SP <- function(betas, sigma, alphas, gama, sigma.b, dat, d., qp, iter, reltol = 1e-08, tol = 1e-06){
                       
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
            mat.bq <- matrix(rep(b[, q], ncol(miss.)), n * (ncol(miss.)), k, TRUE)
            wGH <- as.matrix(expand.grid(lapply(1:q, function(k, u) u$w, u = GH)))
            wGH <- sqrt(2) * apply(wGH, 1, prod) * exp(rowSums(b. * b.))

            sigma.y <- sigma
            gammas <- alphas
            alpha <- gama

            lgLik <- re.warn.miss <- re.warn.miss2 <- re.warn.rand <- numeric(iter)

            environment(opt.missin) <- environment(opt.missin2) <- environment(opt.random) <- environment()
            
    for(i in 1:iter){
        # E-step
        mu.y <- c(X %*% betas) + Z.tbqy
        log.p.yb <- rowsum(dnorm(y, mu.y, sigma.y, log = TRUE), id)
        #log.p.yb <- dnorm(y, mu.y, sigma.y, TRUE) 
        #log.p.yb <- t(sapply(split(log.p.yb, id), function(x) colSums(matrix(x, ncol = k)) ))
        #dimnames(log.p.yb) <- NULL

        mu.t <- plogis(c(W %*% gammas) + alpha * mat.bq)
        mr <- Rn * log(mu.t) + Rn. * log(1 - mu.t)
        prs. <- array(mr, c(ncol(miss.), n, k))
        log.p.tb <- colSums(prs., 2)
    
        p.ytb <- exp(log.p.yb + log.p.tb)
        p.b <- dnorm(b, sd = sigma.b)
        wGHn <- wGH * p.b
        p.yt <- c(p.ytb %*% wGHn)
        p.byt <- p.ytb / p.yt
        post.b <- c(p.byt %*% (bqy * wGHn))
        post.vb <- if(qy == 1) {
            c(p.byt %*% (bqy2 * wGHn)) - c(post.b * post.b)
            } else {
            (p.byt %*% (bqy2 * wGHn)) - t(apply(post.b, 1, function(x) x %o% x))
            }
        lgLik[i] <- sum(log(p.yt))
        cat("iter:", i, "\tlogLik: ", lgLik[i], "\n")
        difL <- lgLik[i] - lgLik[i - 1]
        if(i > 1 && (difL < 0 & abs(difL) > reltol * (abs(lgLik[i - 1]) + reltol)))
            stop("log-likelihood failed to increase.")
        if(i > 1 && difL < tol){
            cat("\nConverged!")
            break
        }
   
        # M-step
        Zb <- if(qy == 1) post.b[id] else rowSums(Z * post.b[id, ])
        betasn <- c(solve(tXX, crossprod(X, y - Zb)))
        mu <- y - c(X %*% betas)
        tr.tZZvarb <- sum(tZZ * post.vb)
        sigman <- sqrt(c(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / N)
    
        thets <- c(gammas)
        optmz.missin <- try(optim(thets, opt.missin, method = "BFGS", control = list(maxit = 5)), TRUE)
        optmz.missin2 <- try(optim(alpha, opt.missin2, method = "BFGS", control = list(maxit = 5)), TRUE)
        
        thets <- log(sigma.b)
        optmz.random <- try(optim(thets, opt.random, method = "BFGS", control = list(maxit = 5)), TRUE)
        
    
        # update parameter values
        betas <- betasn
        sigma.y <- sigman
        if(!inherits(optmz.missin, "try-error")) {
            gammas <- optmz.missin$par
        } else {
            re.warn.miss[i] <- 1
        }    
        if(!inherits(optmz.missin2, "try-error")) {
            alpha <- optmz.missin2$par
        } else {
            re.warn.miss2[i] <- 1
        }    
        if(!inherits(optmz.random, "try-error")) {
            sigma.b <- exp(optmz.random$par)
        } else {
            re.warn.rand[i] <- 1
        }
    }
        warn <- any(c(re.warn.miss, re.warn.miss2, re.warn.rand)==1)
        
        #f.fSP <- function(x){ dnorm(x, 0, sigma.b) }
        #z2 <- apply(grid.., 1, function(x) f.fSP(x))
        #dgrid. <- abs(grid..[2] - grid..[1])
        #ISESP <- sum(dgrid. * (z1 - z2)^2)
        
        if(i==iter) cat("\nNumber of iterations reached without convergence")
        
        list(bts = betas, sig = sigma.y, alps = gammas, gamma. = alpha, 
             vb = (sigma.b)^2, logLik. = lgLik[i], warn = warn, its = i)

}

#SP(betas = cfs.lme, sigma = sigma, alphas = cfs.lmer, gama = 0.6, sigma.b = sigmab, 
#    dat = simul$dat, d. = d., qp = 21, iter = 5000, reltol = 1e-08, tol = 1e-06)

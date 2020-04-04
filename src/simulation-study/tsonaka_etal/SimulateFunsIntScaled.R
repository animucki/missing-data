########################################################################################
# program : SimulateFunsInt.R                                                           
# project : Shared Parameter model with Unknown Mixing distribution                     
# purpose : Functions required to simulate data under the SPM with various univariate   
#           distributions for the random effects (random intercepts only):              
#           1. normal                                                                   
#           2. normal mixture                                                           
#           3. Log-Normal                                                               
#           4. Discrete                                                                 
#                                                                                       
# creation date : 19APRIL2007                                                           
# update :                                                                              
# Note: Data are simulated with zero mean for b                                         
#       Linear mixed model for the measurement process + Mixed effects logistic         
#       regression for the missingness process (Non-monotone missingness)               
#       The 1st measurement is always obseved                                           
# programmer : R. Tsonaka                                                               
########################################################################################


# Simulate SPM with RE from Log-Normal (Random intercepts only)
# sigma is sd
sim.LGSP <- function(n, p, betas, loc, scale., sigma, alphas, gamma., max.time){
                
    # Simulate design matrices
    X <- cbind(1, rep(rbinom(n, size = 1, prob = 0.5), each = p), rep(seq(0, max.time, len = p), n))              
    # Missingness
    W. <- X[, 1:3] # n*p rows
    W <- X[ -seq(from = 1, to = n * p, by = p), 1:3] # n*(p-1) rows we do not allow missingness in the 1st visit
    # Measurement
    X <- cbind(X[, -2], X[, 2] * X[, 3]) # treatment has been excluded as main effect 

    # Simulate random-effects
    b <- rlnorm(n, loc, scale.)
    b <- b - exp(loc + 1/2 * scale.^2)

    # Simulate responses
    # Measurement
    eta.y <- as.vector(X %*% betas + rep(b, each = p))
    y <- rnorm(n * p, eta.y, sigma)
    y <- matrix(y, byrow = TRUE, ncol = p)
                
    # Missingness
    eta.r <- as.vector(W %*% alphas + gamma. * rep(b, each = p - 1))
    pr <- plogis(eta.r)
    r <- rbinom(n * (p - 1), 1, pr)
    r <- matrix(r, byrow = TRUE, ncol = p - 1) == 0

    # Delete missing responses
    y. <- y[, -1]
    y.[r] <- NA
    y[, -1] <- y.
    ind <- !apply(y, 1, function(x) all(is.na(x)))
    y <- y[ind, ]
    #Z <- X[rep(ind, each = p), 1, drop = FALSE]
    #X <- cbind(X[rep(ind, each = p), -1], rep(1:sum(ind), each = p))
    W <- cbind(W[rep(ind, each = p - 1), -1], rep(1:sum(ind), each = p - 1))
    #W <- cbind(W[rep(ind, each = p - 1), -1])
    
list(dat = y, X = cbind(W.[, -1], rep(1:sum(ind), each = p)), W = W, n = n, scale. = scale., loc = loc, b = b, per.mis = sum(is.na(y)) / (n * p))
}
#aa <- sim.LGSP(n = 200, p = 5, betas = c(-1.2, 0.5, -1.5), loc = 0, scale. = 0.97^2, sigma = 0.5, alphas = c(1.6, 2.5, -0.5), gamma. = 0.7, max.time = 5)
#summary(aa$b)
#plot(sort(aa$b), dlnorm(sort(aa$b)))

# Plot of the simulated log-Normal random effects distribution
sim.b.plotLG <- function(prec, sb, loc, n, scale.){
    x <- seq(-sb, sb, length = prec)
    mu. <- exp(loc + 0.5 * scale.^2)
    grid. <- as.matrix(x)
    f <- function(x){ 
        dens <- dlnorm(x, loc, scale.)
        }    
    z1 <- apply(grid., 1, function(x) f(x))
    dim(z1) <- c(prec)
    #par(mfrow = c(2,2))
    plot(x - mu., z1, lwd = 2, col = 1, cex.lab = 2, cex.axis = 1.9, type = "l", xlab = "", ylab = "", xlim = c(-10, 10), ylim = c(0, 0.7))
    mtext("Random Intercepts", side = 1, line = 4, cex = 1.9)
    mtext("Density", side = 2, line = 4, cex = 1.9)
}
#sim.b.plotLG(prec = 100, sb = 10, loc = 0, n = 100, scale. = 0.97^2)


# Simulate SPM with a Normal RE (random intercepts)
# sigmab is sd
sim.normal <- function(n, p, betas, meanb, sigmab, sigma, alphas, gamma., max.time){
                     
    # Simulate design matrices
    X <- cbind(1, rep(rbinom(n, size = 1, prob = 0.5), each = p), rep(seq(0, max.time, length = p), n))
    # Missingness
    W. <- X[, 1:3]
    W <- X[ -seq(from = 1, to = n * p, by = p), 1:3]
    # Measurement
    X <- cbind(X[, -2], X[, 2] * X[, 3]) # treatment has been excluded as main effect 
                
    # Simulate random-effects
    b <- rnorm(n, 0, sigmab)

    # Simulate responses
    # Measurement
    eta.y <- as.vector(X %*% betas + rep(b, each = p))
    y <- rnorm(n * p, eta.y, sigma)
    y <- matrix(y, byrow = TRUE, ncol = p)
                
    # Missingness
    eta.r <- as.vector(W %*% alphas + gamma. * rep(b, each = p - 1))
    pr <- plogis(eta.r)
    r <- rbinom(n * (p - 1), 1, pr)
    r <- matrix(r, byrow = TRUE, ncol = p - 1) == 0

    # Delete missing responses
    y. <- y[, -1]
    y.[r] <- NA
    y[, -1] <- y.
    ind <- !apply(y, 1, function(x) all(is.na(x)))
    y <- y[ind, ]
    #Z <- X[rep(ind, each = p), 1, drop = FALSE]
    #X <- cbind(X[rep(ind, each = p), -1], rep(1:sum(ind), each = p))
    W <- cbind(W[rep(ind, each = p - 1), -1], rep(1:sum(ind), each = p - 1))
    #W <- cbind(W[rep(ind, each = p - 1), -1])

list(dat = y, X = cbind(W.[, -1], rep(1:sum(ind), each = p)), W = W, n = n, sigmab = sigmab, meanb = meanb, b = b, per.mis = sum(is.na(y)) / (n * p))
}
#aa <- sim.normal(n = 200, p = 5, betas =  c(-1.2, 0.5, -1.5), meanb = 0, sigmab = 2, sigma = 0.5, alphas = c(1.6, 2.5, -0.5), gamma. = 0.7, max.time = 5)


# Plot of the simulated normal random effects distribution SPM
sim.b.plotN <- function(prec, sb, n, D.Y){
    x <- seq(-sb, sb, length = prec)
    grid. <- as.matrix(x)
    f <- function(x){
        mu. <- 0
        dens <- dnorm(x, mu., D.Y)
    }
    z1 <- apply(grid., 1, function(x) f(x))
    dim(z1) <- c(prec)
    plot(x, z1, col = 1, cex = 0.45, type = "l", lwd = 2)
}
#sim.b.plotN(prec = 100, sb = 10, n = 200, D.Y = 2)


# Simulate SPM with a mixture of normals RE (random intercepts)
sim.mix <- function(n, p, betas.c, n.c, probs, betas, meanb, sigmab, sigma, alphas, gamma., max.time){
                     
    # Simulate design matrices
    X <- cbind(1, rep(rbinom(n, size = 1, prob = 0.5), each = p), rep(seq(0, max.time, length = p), n))
    # Missingness
    W. <- X[, 1:3]
    W <- X[ -seq(from = 1, to = n * p, by = p), 1:3]
    # Measurement
    X <- cbind(X[, -2], X[, 2] * X[, 3]) # treatment has been excluded as main effect 
                
    # Simulate random-effects
    ind <- sample(1:n.c, n, TRUE, probs)
    n.ind <- as.vector(table(ind))
    b <- numeric(n)
    for(i in 1:n.c){
        b[ind == i] <- rnorm(n.ind[i], betas.c[i], sigmab)
    }
    b <- b - sum(betas.c * probs)

    # Simulate responses
    # Measurement
    eta.y <- as.vector(X %*% betas + rep(b, each = p))
    y <- rnorm(n * p, eta.y, sigma)
    y <- matrix(y, byrow = TRUE, ncol = p)
                
    # Missingness
    eta.r <- as.vector(W %*% alphas + gamma. * rep(b, each = p - 1))
    pr <- plogis(eta.r)
    r <- rbinom(n * (p - 1), 1, pr)
    r <- matrix(r, byrow = TRUE, ncol = p - 1) == 0

    # Delete missing responses
    y. <- y[, -1]
    y.[r] <- NA
    y[, -1] <- y.
    ind <- !apply(y, 1, function(x) all(is.na(x)))
    y <- y[ind, ]
    #Z <- X[rep(ind, each = p), 1, drop = FALSE]
    #X <- cbind(X[rep(ind, each = p), -1], rep(1:sum(ind), each = p))
    W <- cbind(W[rep(ind, each = p - 1), -1], rep(1:sum(ind), each = p - 1))
    #W <- cbind(W[rep(ind, each = p - 1), -1])

list(dat = y, X = cbind(W.[, -1], rep(1:sum(ind), each = p)), W = W, n = n, sigmab = sigmab, meanb = meanb, b = b, per.mis = sum(is.na(y)) / (n * p))
}


# Plot of the simulated log-Normal random effects distribution
sim.b.plotMIX <- function(prec, sb, betas.c, n.c, probs, n, sigmab){
    x <- seq(-sb, sb, length = prec)
    mu. <- sum(betas.c * probs)
    grid. <- as.matrix(x)
    f <- function(x){ 
        dens <- sum(dnorm(x, betas.c, sigmab) * probs)
        }    
    z1 <- apply(grid., 1, function(x) f(x))
    dim(z1) <- c(prec)
    #par(mfrow = c(2,2))
    plot(x - mu., z1, col = 1, lwd = 2, cex.lab = 2, cex.axis = 1.9, type = "l", xlab = "", ylab = "", xlim = c(-10, 10), ylim = c(0, 0.7))
    #mtext("Random Intercepts", side = 1, line = 4, cex = 1.9)
    #mtext("Density", side = 2, line = 4, cex = 1.9)
}
#sim.b.plotMIX(prec = 200, sb = 10, betas.c =c(-2.2, -0.05, 0.05, 2.2), n.c = 4, probs = c(0.4, 0.1, 0.1, 0.4), n = 100, sigmab = 0.4)


# Simulate SPM with a discrete RE (random intercepts)
sim.discr <- function(n, p, betas.c, n.c, probs, betas, meanb, sigmab, sigma, alphas, gamma., max.time){

    # Simulate design matrices
    X <- cbind(1, rep(rbinom(n, size = 1, prob = 0.5), each = p), rep(seq(0, max.time, length = p), n))
    # Missingness
    W. <- X[, 1:3]
    W <- X[ -seq(from = 1, to = n * p, by = p), 1:3]
    # Measurement
    X <- cbind(X[, -2], X[, 2] * X[, 3]) # treatment has been excluded as main effect 
                
    # Simulate random-effects
    ind <- sample(1:n.c, n, TRUE, probs)
    n.ind <- as.vector(table(ind))
    b <- numeric(n)
    for(i in 1:n.c){
        b[ind == i] <- rep(betas.c[i], n.ind[i])
    }
    b <- b - sum(betas.c * probs)

    # Simulate responses
    # Measurement
    eta.y <- as.vector(X %*% betas + rep(b, each = p))
    y <- rnorm(n * p, eta.y, sigma)
    y <- matrix(y, byrow = TRUE, ncol = p)
                
    # Missingness
    eta.r <- as.vector(W %*% alphas + gamma. * rep(b, each = p - 1))
    pr <- plogis(eta.r)
    r <- rbinom(n * (p - 1), 1, pr)
    r <- matrix(r, byrow = TRUE, ncol = p - 1) == 0

    # Delete missing responses
    y. <- y[, -1]
    y.[r] <- NA
    y[, -1] <- y.
    ind <- !apply(y, 1, function(x) all(is.na(x)))
    y <- y[ind, ]
    #Z <- X[rep(ind, each = p), 1, drop = FALSE]
    #X <- cbind(X[rep(ind, each = p), -1], rep(1:sum(ind), each = p))
    W <- cbind(W[rep(ind, each = p - 1), -1], rep(1:sum(ind), each = p - 1))
    #W <- cbind(W[rep(ind, each = p - 1), -1])

list(dat = y, X = cbind(W.[, -1], rep(1:sum(ind), each = p)), W = W, n = n, sigmab = sigmab, meanb = meanb, b = b, per.mis = sum(is.na(y)) / (n * p))
}
#aa <- sim.normal(n = 100, p = 5, betas =  c(-1.2, 0.5, 1.5), meanb = 0, sigmab = 2, sigma = 0.5, alphas = c(1.6, 2.5, -0.5), gamma. = 0.7)


# Plot of the simulated discrete random effects distribution
sim.b.plotDISCR <- function(prec, sb, betas.c, n.c, probs, n){
    plot(betas.c, probs, col = 1, cex.lab = 2, cex.axis = 1.9, type = "n", xlab = "", ylab = "", xlim = c(-sb, sb))
    segments(x0 = betas.c, y0 = rep(0, length(betas.c)), 
             x1 = betas.c, y1 = probs, lwd = 2)

}
#sim.b.plotDISCR(prec = 100, sb = 6, betas.c = c(-2.5, -1.5, 1.5, 2.5), n.c =4, probs = c(0.3, 0.2, 0.2, 0.3), n = 100)


# Plot of the fitted random effects distribution - NPSP
fitted.b.plot <- function(prec, sb, sup, D.Y, typ = "n"){
            
    plot(sup[, 1], sup[, 2], col = 1, cex.lab = 2, cex.axis = 1.9, type = typ, xlab = "", ylab = "", xlim = c(-sb, sb))
    segments(x0 = sup[, 1], y0 = rep(0, length(sup[, 1])), 
             x1 = sup[, 1], y1 = sup[, 2])
    mtext("Grid Points", side = 1, line = 4, cex = 1.9)
    mtext("Weights", side = 2, line = 4, cex = 1.9)

}

#fitted.b.plot(prec = 10, sb = 5, sup = model.1$sup.points, D.Y = 1, typ = "n")

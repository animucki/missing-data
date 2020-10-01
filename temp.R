
mu2 <- c(-1, 0.2)
sigma2 <- matrix(c(1.6, 1.2, 1.2, 1.4), ncol = 2)
library(mc2d)

a <- dlm::arms(y.start = c(0,0),
               myldens = function (x) -(x %*% solve(sigma2) %*% x),
               indFunc = function (x) {sum(abs(x)) < 10},
               n.sample = 100)
plot(a)

b <- dlm::arms(y.start = 0,
               myldens = function (x) - x - x^2,
               indFunc = function (x) {x >= -5 & x < 5},
               n.sample = 1e6)
#hist(b)
plot(density(b, from = -5, to = 5))
mean(b)


#Copula
library(mc2d)

n <- 1000
rho <- 0.95
sigmab <- 2.5

corMat <- matrix(c(1,rho,rho,1),ncol = 2)
#aMat <- chol(corMat)
#zMat <- rmultinormal(n, c(0,0))
#xMat <- zMat %*% aMat
xMat <- rmultinormal(n, sigma = c(1,rho,rho,1))
plot(xMat)

uMat <- pnorm(xMat)
plot(uMat)
hist(uMat[,1])
hist(uMat[,2])

bMat <- 2 * sqrt(3) * sigmab * (uMat - 0.5)
plot(bMat)
sqrt(diag(var(bMat)))
cor(bMat)


erf <- function (x) {
  2*(pnorm(x) - 0.5)
}

erfPrime <- function (x) {
  2 * dnorm(x, sd = sqrt(2))
}
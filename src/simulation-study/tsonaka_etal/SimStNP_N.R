###########################################################################################
# program : SimStNP_N.R                                                            
# project : Shared Parameter models with Flexible Mixing distributions (Non Parametric)    
# purpose : Simulation study to compare the NPSP model fitted under VEM                    
#           with the simple SP model                                                       
# require: SupportFuns.R                                                                   
# creation date : 06DEC2007                                                                
# programmer : R. Tsonaka                                                                  
###########################################################################################

#########################################
#                                       #
#  Simulate from Normal                 #
#                                       #
#########################################


source(paste(file., "SimulateFunsIntScaled.R", sep = ""))
source(paste(file., "SupportFunsIntScaled.R", sep = ""))

library(nlme)
library(lme4)

# initial values
in.pars <- list(n = 200, p = 5, betas = c(-1.2, 0.5, -1.5), meanb = 0, sigma = sqrt(0.5), 
                sigmab = sqrt(2), alphas = c(1.6, 2.5, -0.5), gamma. = 1.5)

set.seed(123)

warnSPSP <- 0
warnSP <- 0
      
    simul <- sim.normal(n = in.pars$n, p = in.pars$p, betas = in.pars$betas, meanb = in.pars$meanb, 
                      sigmab = in.pars$sigmab, sigma = in.pars$sigma, alphas = in.pars$alphas, 
                      gamma. = in.pars$gamma., max.time = in.pars$p) # sigma is stdev
    summary(simul$b)
    miss <- t(apply(simul$dat, 1, function(x){ x <- ifelse(!is.na(x), 1, 0); x}))
    1-sum(miss)/1000
    d. <- data.frame(response = c(t(simul$dat)), miss = c(t(miss)), 
                     treat = simul$X[, 1], time. = simul$X[, 2], id = simul$X[, 3])
    d.lmer <- data.frame(miss = c(t(miss[, -1])), treat = simul$W[, 1], time. = simul$W[, 2], id = simul$W[, 3])
    model.lme <- lme(response ~ time. * treat - treat, random = ~ 1|id, data = d., na.action = na.exclude)                 
    cfs.lme <- model.lme$coef$fixed
    cfs.lme
    names(cfs.lme) <- NULL
    sigma <- model.lme$sigma
    sigma
    sigmab <- sqrt(unlist(lapply(pdMatrix(model.lme$modelStruct$reStruct), "*", model.lme$sigma^2)))
    sigmab
    X <- model.matrix(~ time. * treat - treat - 1, data = d.)
    attr(X, "contrasts") <- NULL
    attr(X, "assign") <- NULL
    model.lmer <- glmer(miss ~ treat + time. + (1 | id), data = d.lmer, 
                       nAGQ=1, family = binomial(link = "logit"), na.action = na.exclude)  
    W <- model.matrix(~ treat + time., data = d.lmer)
    attr(W, "contrasts") <- NULL
    attr(W, "assign") <- NULL
    cfs.lmer <- fixef(model.lmer)
    cfs.lmer
    simul$per 
    miss <- miss[, -1]

NPst.time <- proc.time()      
    
    # SPSP model
    cat("\nModel SPSP")
    SPSP.error <- try({        
    model.1 <- NPSP(data. = simul$dat, miss. = miss, X = X, W = W, betas = in.pars$betas, sigma2 = in.pars$sigma^2, 
                    sigmab = sigmab, alphas = in.pars$alphas, delta = in.pars$gamma., gp = 161, bound = 4, 
                    reltol = 1e-08, epsilon = 0.001, iter = 2000)
                }, TRUE)

    if(!inherits(SPSP.error, "try-error")){
        RES.SPSP <- model.1} else {
        warnSPSP <- 1
        }  

NPtot.time <- proc.time() - NPst.time


SPst.time <- proc.time()                  
         
    # SP model
    cat("\nModel SP")
    SP.error <- try({
    model.2 <- SP(betas = in.pars$betas, sigma = in.pars$sigma, alphas = in.pars$alphas, gama = in.pars$gamma., sigma.b = in.pars$sigmab, 
                  dat = simul$dat, d. = d., qp = 21, iter = 2000, reltol = 1e-08, tol = 1e-07)
                }, TRUE)

        if(!inherits(SP.error, "try-error")){
                    RES.SP <- model.2} else {
                    warnSP <- 1
                    }        
   
SPtot.time <- proc.time() - SPst.time

NPbetas <- model.1$bts
NPalphas <- model.1$alps
NPalphas2 <- model.1$alphas
NPgamma <- model.1$gamma.
NPvb <- model.1$vb
NPmb <- model.1$mb
NPsig <- sqrt(model.1$sig)
NPmiss <- model.1$per.mis
NPsb <- model.1$sigmab
               
SPbetas <- model.2$bts
SPalphas <- model.2$alps
SPgamma <- model.2$gamma.
SPvb <- model.2$vb
SPsig <- model.2$sig


# Std Errors
source(paste(file.., "SupportFunsStdErrorsScaled.R", sep = ""))


###############
# SPSP
Sup <- model.1$sup.points
pi. <- Sup[, 2]
gridp <- Sup[, 1] - NPmb

# Profile
data. <- simul$dat
y <- c(t(data.))
n <- nrow(data.) # sample size
p <- ncol(data.)
na <- ncol(W)
miss. <- miss
gp <- length(gridp)

hs <- hess.parsProfile(c(c(NPbetas, NPsig^2, NPalphas, NPgamma)))
file.out <- paste(file., "NPHessianProfile", jobnr, ".txt", sep = "") 
write.table(hs, file.out)   


gridp <- Sup[, 1] - NPmb 
pis <- Sup[, 2]
gp <- length(gridp)

hs <- hess.parsThres(c(c(NPbetas, NPsig^2, NPalphas, NPgamma, pis[-1])))
file.out <- paste(file., "NPHessian", jobnr, ".txt", sep = "") 
write.table(hs, file.out)   



###############
# SPM
dat <- data.
qp <- 21
hs <- hess.pars(c(c(SPbetas, SPsig, SPalphas, SPgamma, sqrt(SPvb))))
file.out <- paste(file., "SPHessian", jobnr, ".txt", sep = "") 
write.table(hs, file.out)   

# Grouped 0.05
thres <- 0.05
indh <- which(Sup[, 2] < thres)
gridp <- Sup[-indh, 1] - NPmb 
pis <- Sup[-indh, 2]
pr.pis <- pis / sum(pis)
pis.new <- pis + (1 - sum(pis)) * pr.pis
gp <- length(gridp)
hs <- hess.parsThres(c(c(NPbetas, NPsig^2, NPalphas, NPgamma, pis.new[-1])))
file.out <- paste(file., "NPHessian005Thres", jobnr, ".txt", sep = "") 
write.table(hs, file.out)   

# Grouped
thres <- 0.01
indh <- which(Sup[, 2] < thres)
gridp <- Sup[-indh, 1] - NPmb 
pis <- Sup[-indh, 2]
pr.pis <- pis / sum(pis)
pis.new <- pis + (1 - sum(pis)) * pr.pis
gp <- length(gridp)
hs <- hess.parsThres(c(c(NPbetas, NPsig^2, NPalphas, NPgamma, pis.new[-1])))
file.out <- paste(file., "NPHessian001Thres", jobnr, ".txt", sep = "") 
write.table(hs, file.out)   

# Grouped
thres <- 0.02
indh <- which(Sup[, 2] < thres)
gridp <- Sup[-indh, 1] - NPmb 
pis <- Sup[-indh, 2]
pr.pis <- pis / sum(pis)
pis.new <- pis + (1 - sum(pis)) * pr.pis
gp <- length(gridp)
hs <- hess.parsThres(c(c(NPbetas, NPsig^2, NPalphas, NPgamma, pis.new[-1])))
file.out <- paste(file., "NPHessian002Thres", jobnr, ".txt", sep = "") 
write.table(hs, file.out)   




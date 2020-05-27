# clear workspace
rm(list = ls())

library(armspp)
library(futile.logger)
library(lme4)
library(mc2d)
library(numDeriv)
library(tictoc)
library(tidyverse)

# flog.appender(appender.file(paste0('./log/emp-', format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss"), ".log"))) %>% invisible
flog.threshold('trace') %>% invisible

flog.info('Sourcing functions...')
source("src/common/utils.r")
source("src/empirical-analysis/fit.parametric.r")
source("src/empirical-analysis/fit.hybrid.r")
source("src/empirical-analysis/fit.class.r")
source("src/empirical-analysis/fit.npsp.r")


flog.info('Loading data...')
dataWide <- read.csv('./data/scleroderma.csv') %>% filter(is.finite(MAXFIB))

#Make data into long format
data <- dataWide %>%
  pivot_longer(cols = 4:11,
               names_to = "time",
               names_prefix = "FVC",
               values_to = "FVC") %>%
  mutate(time = as.integer(time),
         CYC = as.integer(TXGROUP == 'A'),
         observed = as.integer(!is.na(FVC)))

#Ignorable analysis
if(file.exists("data/m1.rds") && file.exists("data/m1r.rds")) {
  flog.info("Reading ignorable model from file...")
  m1 <- readRDS("data/m1.rds")
  m1r <- readRDS("data/m1r.rds")
} else {
  flog.info("Fitting ignorable model...")
  m1 <- lmer(FVC ~ (1 | pt_id) + CYC * (time + FVC0 + MAXFIB), data = data, REML = F)
  m1r <- glm(observed ~  CYC + time + FVC0 + MAXFIB, data = data, family = binomial())
  #m1r <- glm(observed ~  1, data = data, family = binomial())
  saveRDS(m1, "data/m1.rds")
  saveRDS(m1r, "data/m1r.rds")
}

y <- dataWide %>% select(starts_with("FVC")) %>% select(-1) %>% as.matrix
r <- as.integer(!is.na(y))
dim(r) <- dim(y)

X <- data %>% transmute("(Intercept)"=1, CYC = as.integer(TXGROUP == 'A'), time = as.integer(time), FVC0 = FVC0, MAXFIB = MAXFIB,
                        "CYC:time" = CYC*time, "CYC:FVC0" = CYC*FVC0, "CYC:MAXFIB" = CYC*MAXFIB) %>% as.matrix

W <- data %>% transmute("(Intercept)"=1, CYC = as.integer(TXGROUP == 'A'), time = as.integer(time), FVC0 = FVC0, MAXFIB = MAXFIB) %>% as.matrix

init <- c(fixef(m1), as.data.frame(VarCorr(m1))$sdcor[1], sigma(m1))

#SPM
if(file.exists("data/m2.rds")) {
  flog.info("Reading shared-parameter model from file...")
  m2 <- readRDS("data/m2.rds")
} else {
  flog.info("Fitting shared-parameter model...")
  set.seed(1234L)
  m2 <- fit.parametric(y = y, r = r, X = X, W = W, init = init)
  saveRDS(m2, "data/m2.rds")
}
flog.info("SPM successfully loaded!")

# Hybrid
if(file.exists("data/m3.rds")) {
  flog.info("Reading hybrid model from file...")
  m3 <- readRDS("data/m3.rds")
} else {
  flog.info("Fitting hybrid model...")
  set.seed(111L)
  m3 <- fit.hybrid(y=y, r=r, X=X, W=W, nClasses=3, init=init)
  saveRDS(m3, "data/m3.rds")
}
flog.info("Hybrid model successfully loaded!")

# Class/survival
if(file.exists("data/m4.rds")) {
  flog.info("Reading class/survival model from file...")
  m4 <- readRDS("data/m4.rds")
} else {
  flog.info("Fitting class/survival model...")
  set.seed(1121L)

  isExcl <- rowSums(r) <= 1
  rowsExcl <- rep(isExcl, each = ncol(r))

  m4 <- fit.class(y=y, r=r, X=X, W=W[,-1], nClasses=3, init=init)

  saveRDS(m4, "data/m4.rds")
}
flog.info("Class/survival model successfully loaded!")

# SPSPM (Tsonaka et. al)
if(file.exists("data/m5.rds")) {
  flog.info("Reading SPSP model from file...")
  m5 <- readRDS("data/m5.rds")
} else {
  flog.info("Fitting SPSP model...")
  set.seed(1151L)
  m5 <- fit.npsp(y=y, r=r, X=X, W=W, init=init)
  saveRDS(m5, "data/m5.rds")
}
flog.info("Semi-parametric SPM successfully loaded!")

# ----- Fitting null models -----
#Ignorable analysis
if(file.exists("data/m1_null.rds")) {
  flog.info("Reading (null) ignorable model from file...")
  m1_null <- readRDS("data/m1_null.rds")
} else {
  flog.info("Fitting (null) ignorable model...")
  m1_null <- lmer(FVC ~ (1 | pt_id) + time + FVC0 + MAXFIB, data = data, REML = F)
  saveRDS(m1_null, "data/m1_null.rds")
}

X_null <- data %>% transmute("(Intercept)"=1, time = as.integer(time), FVC0 = FVC0, MAXFIB = MAXFIB) %>% as.matrix
init_null <- c(fixef(m1_null), as.data.frame(VarCorr(m1_null))$sdcor[1], sigma(m1_null))

#SPM
if(file.exists("data/m2_null.rds")) {
  flog.info("Reading shared-parameter model from file...")
  m2_null <- readRDS("data/m2_null.rds")
} else {
  flog.info("Fitting shared-parameter model...")
  set.seed(12384L)
  m2_null <- fit.parametric(y = y, r = r, X = X_null, W = W, init = init_null)
  saveRDS(m2_null, "data/m2_null.rds")
}
flog.info("SPM (null) successfully loaded!")

# Hybrid
if(file.exists("data/m3_null.rds")) {
  flog.info("Reading hybrid model from file...")
  m3_null <- readRDS("data/m3_null.rds")
} else {
  flog.info("Fitting hybrid model...")
  set.seed(1118L)
  m3_null <- fit.hybrid(y=y, r=r, X=X_null, W=W, nClasses=3, init=init_null)
  saveRDS(m3_null, "data/m3_null.rds")
}
flog.info("Hybrid model (null) successfully loaded!")

# Class/survival
if(file.exists("data/m4_null.rds")) {
  flog.info("Reading class/survival model from file...")
  m4_null <- readRDS("data/m4_null.rds")
} else {
  flog.info("Fitting class/survival model...")
  set.seed(11821L)

  m4_null <- fit.class(y=y, r=r, X=X_null, W=W[,-1], nClasses=3, init=init_null)

  saveRDS(m4_null, "data/m4_null.rds")
}
flog.info("Class/survival model (null) successfully loaded!")

# SPSPM (Tsonaka et. al)
if(file.exists("data/m5_null.rds")) {
  flog.info("Reading SPSP model from file...")
  m5_null <- readRDS("data/m5_null.rds")
} else {
  flog.info("Fitting SPSP model...")
  set.seed(11851L)
  m5_null <- fit.npsp(y=y, r=r, X=X_null, W=W, init=init_null)
  saveRDS(m5_null, "data/m5_null.rds")
}
flog.info("Semi-parametric SPM (null) successfully loaded!")

#---- create comparison tables

models <- list(m1, m2, m3, m4, m5)
models_null <- list(m1_null, m2_null, m3_null, m4_null, m5_null)
modelNames <- c("Ignorable", "SPM", "Hybrid", "Class", "SPSP")

#---- first: likelihood table
table1 <- lapply(models[-1],
                 function(m) data.frame(neg2ll = m$res,
                                        aic = m$res + 2*length(unlist(m$pars)),
                                        bic = m$res + log(nrow(y)) * length(unlist(m$pars))
                 )) %>% bind_rows
table1 <- bind_rows(data.frame(neg2ll = deviance(m1) + deviance(m1r)) %>%
                      mutate(aic = neg2ll + 2 * (10+5),
                             bic = neg2ll + log(nrow(y)) * (10+5)),
                    table1)
rownames(table1) <- modelNames
table1

#---- second: parameters table
table2 <- list()
for (i in seq_along(models)) {
  if(i==1) {
    pars <- fixef(m1)[-1]
    se <- sqrt(diag(vcov(m1)))[-1]
    names(se) <- paste0("se_", names(pars))
  } else if(i==5) {
    pars <- m5$pars[2:8]
    se <- sqrt(diag(solve(m5$hess)))[2:8]

    names(pars) <- colnames(X)[2:8]
    names(se) <- paste0("se_", names(pars))

  } else {
    pars <- models[[i]]$pars$beta[2:8]
    se <- sqrt(diag(solve(models[[i]]$hess)))[2:8]

    names(pars) <- colnames(X)[2:8]
    names(se) <- paste0("se_", names(pars))
  }

  table2[[i]] <- as.list(c(pars, se, m=modelNames[i]))
}
table2 <- bind_rows(table2)

order <- c("time", "FVC0", "MAXFIB", "CYC", "CYC:FVC0", "CYC:MAXFIB", "CYC:time")

t2a <- as.numeric(t(table2[,order]))
dim(t2a) <- dim(t(table2[,order]))
t2b <- as.numeric(t(table2[,paste0("se_",order)]))

signif2 <- abs(t2a/t2b) > qnorm(0.975)
signif2

t2a <- format(t2a, digits = 2)
colnames(t2a) <- modelNames
rownames(t2a) <- order

t2b <- paste0("(",format(t2b, digits = 2),")")
dim(t2b) <- dim(t(table2[,paste0("se_",order)]))
colnames(t2b) <- modelNames
rownames(t2b) <- paste0("se_",order)

orderInterleaved <- as.vector(matrix(c(order, paste0("se_",order)), nrow = 2, byrow = T))
table2 <- rbind(t2a, t2b)[orderInterleaved,]
table2

#---- third: "overall effect at time" table
timesT <- c(6,9,12,15,18)
table3se <- table3est <- matrix(NA_real_, nrow = length(timesT), ncol = length(models))

Xt <- as.data.frame(X[X[,"time"] %in% timesT, c("(Intercept)","time","FVC0","MAXFIB")])
Xt0 <- Xt %>% mutate(CYC = 0, cyctime = 0, cycfvc0 = 0, cycmaxfib = 0) %>% as.matrix
Xt1 <- Xt %>% mutate(CYC = 1, cyctime = time, cycfvc0 = FVC0, cycmaxfib = MAXFIB) %>% as.matrix

for (m in seq_along(models)) {
  if (m==1) {
    y0 <- Xt0 %*% fixef(m1)[c("(Intercept)","time","FVC0","MAXFIB", "CYC", "CYC:time", "CYC:FVC0", "CYC:MAXFIB")]
    y1 <- Xt1 %*% fixef(m1)[c("(Intercept)","time","FVC0","MAXFIB", "CYC", "CYC:time", "CYC:FVC0", "CYC:MAXFIB")]
  } else if (m==2) {
    y0 <- Xt0 %*% m2$pars$beta[c(1,3,4,5,2,6,7,8)]
    y1 <- Xt1 %*% m2$pars$beta[c(1,3,4,5,2,6,7,8)]
  } else if (m==3 || m==4) {
    y0 <- y1 <- rep(0, nrow(Xt))
    #expectation over groups (since the group intercept does not have zero expectation)
    for (k in 1:3) {
      y0 <- y0 + ( Xt0 %*% models[[m]]$pars$beta[c(1,3,4,5,2,6,7,8)] + c(0, models[[m]]$pars$mu)[k] ) * softmax(c(0, models[[m]]$pars$eta))[k]
      y1 <- y1 + ( Xt1 %*% models[[m]]$pars$beta[c(1,3,4,5,2,6,7,8)] + c(0, models[[m]]$pars$mu)[k] ) * softmax(c(0, models[[m]]$pars$eta))[k]
    }
  } else if (m==5) {
    y0 <- Xt0 %*% m5$pars[paste0("beta.",c("(Intercept)","time","FVC0","MAXFIB", "CYC", "CYC:time", "CYC:FVC0", "CYC:MAXFIB"))]
    y1 <- Xt1 %*% m5$pars[paste0("beta.",c("(Intercept)","time","FVC0","MAXFIB", "CYC", "CYC:time", "CYC:FVC0", "CYC:MAXFIB"))]
  }

  for (t in seq_along(timesT)) {
    table3est[t, m] <- mean( (y1 - y0)[Xt[,"time"]==timesT[t]] )
    table3se[t, m] <- sd( (y1 - y0)[Xt[,"time"]==timesT[t]] )/sqrt(nrow(y))
  }
}


t3a <- format(table3est, digits = 2)
rownames(t3a) <- paste(timesT,"months")
t3b <- paste0("(",format(table3se, digits = 2),")")
dim(t3b) <- dim(table3se)
table3 <- rbind(t3a, t3b)[as.vector(matrix(1:10, ncol = 5, byrow = T)),]
colnames(table3) <- modelNames
table3

#significances at alpha=5%
abs(table3est/table3se) > qnorm(0.975)

#---- fourth: LRT p-values table
table4 <- list()
for (i in seq_along(models)) {
  if(i==1) {
    table4[[i]] <- data.frame( stat = deviance(m1_null) - deviance(m1), p = pchisq(deviance(m1_null) - deviance(m1), df = 4, lower.tail = F) )
  } else {
    table4[[i]] <- data.frame( stat = models_null[[i]]$res - models[[i]]$res, p = pchisq(models_null[[i]]$res - models[[i]]$res, df = 4, lower.tail = F) )
  }
}
table4 <- bind_rows(table4)
rownames(table4) <- modelNames
format(table4, nsmall=4, digits=0, scientific = F)
# p-values obtained via likelihood ratio test.
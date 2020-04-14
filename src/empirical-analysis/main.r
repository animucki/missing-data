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
# source("src/empirical-analysis/fit.class.r")
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
         CYC = as.integer(TXGROUP == 'A'))

#Ignorable analysis
if(file.exists("data/m1.rds")) {
  flog.info("Reading ignorable model from file...")
  m1 <- readRDS("data/m1.rds")
} else {
  flog.info("Fitting ignorable model...")
  m1 <- lmer(FVC ~ (1 | pt_id) + CYC * (time + FVC0 + MAXFIB), data = data, REML = F)
  saveRDS(m1, "data/m1.rds")
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
# if(file.exists("data/m4.rds")) {
#   flog.info("Reading class/survival model from file...")
#   m4 <- readRDS("data/m4.rds")
# } else {
#   flog.info("Fitting class/survival model...")
#   set.seed(1121L)
#   m4 <- fit.class(y=y, r=r, X=X, W=W[,-1], nClasses=3, init=init)
#   saveRDS(m4, "data/m4.rds")
# }
# flog.info("Class/survival model successfully loaded!")

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
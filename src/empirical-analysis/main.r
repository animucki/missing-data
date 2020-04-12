# clear workspace
rm(list = ls())

library(armspp)
library(futile.logger)
library(lme4)
library(mc2d)
library(numDeriv)
library(optimx)
library(tictoc)
library(tidyverse)

# flog.appender(appender.file(paste0('./log/emp-', format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss"), ".log"))) %>% invisible
flog.threshold('trace') %>% invisible

flog.info('Sourcing functions...')
source("src/common/utils.r")
source("src/empirical-analysis/fit.parametric.r")
source("src/empirical-analysis/fit.hybrid.r")


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
m1 <- lmer(FVC ~ (1|pt_id) + CYC*(time + FVC0 + MAXFIB), data = data, REML = F)

y <- dataWide %>% select(starts_with("FVC")) %>% select(-1) %>% as.matrix
r <- as.integer(!is.na(y))
dim(r) <- dim(y)

X <- data %>% transmute("(Intercept)"=1, CYC = as.integer(TXGROUP == 'A'), time = as.integer(time), FVC0 = FVC0, MAXFIB = MAXFIB,
                        "CYC:time" = CYC*time, "CYC:FVC0" = CYC*FVC0, "CYC:MAXFIB" = CYC*MAXFIB) %>% as.matrix

W <- data %>% transmute("(Intercept)"=1, CYC = as.integer(TXGROUP == 'A'), time = as.integer(time), FVC0 = FVC0, MAXFIB = MAXFIB) %>% as.matrix

init <- c(fixef(m1), as.data.frame(VarCorr(m1))$sdcor[1], sigma(m1))

#SPM
set.seed(1234L)
# m2 <- fit.parametric(y=y, r=r, X=X, W=W, init=init)

set.seed(111L)
m3 <- fit.hybrid(y=y, r=r, X=X, W=W, nClasses=3, init=init)
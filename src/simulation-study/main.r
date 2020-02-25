# clear workspace
rm(list = ls())

library(fastGHQuad)
library(futile.logger)
library(lme4)
library(mc2d)
library(mvtnorm)
library(numDeriv)
library(parallel)
library(tictoc)
library(tidyverse)

options(mc.cores = parallel::detectCores() - 1)

# flog.appender(appender.file(paste0('./log/', format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss"), ".log")))
flog.threshold('trace')

flog.info('Sourcing functions...')
source("src/simulation-study/utils.r")
source("src/simulation-study/generateSamples.r")
source("src/simulation-study/fit.ignorable.r")

source("src/simulation-study/fit.parametric.r")
source("src/simulation-study/fit.parametric2.r")

source("src/simulation-study/fit.hybrid.r")
source("src/simulation-study/fit.hybrid2.r")

source("src/simulation-study/fit.class.r")
source("src/simulation-study/fit.class2.r")

set.seed(666L)
df1 <- generateSamples(samples = 1, participants = 100)

# ALL SAMPLES
res <- list()

flog.info('Testing new parametric model...')
tic()
# res[[1]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% lapply(fit.parametric2) %>% bind_rows
res[[1]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% lapply(fit.hybrid2) %>% bind_rows
toc()

# flog.info('Fitting models to MAR scenario...')
# res[[1]] <- df1 %>% mutate(y=yMAR, r=rMAR) %>% group_split(sample) %>% mclapply(fit.multiple) %>% bind_rows %>% mutate(scenario='MAR')

# flog.info('Fitting models to MNAR scenario...')
# res[[2]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% mclapply(fit.multiple) %>% bind_rows %>% mutate(scenario='MNAR')

result <- bind_rows(res)


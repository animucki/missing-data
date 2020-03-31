# clear workspace
rm(list = ls())

library(armspp)
library(fastGHQuad)
library(futile.logger)
library(lme4)
library(mc2d)
library(numDeriv)
library(parallel)
library(tictoc)
library(tidyverse)

options(mc.cores = parallel::detectCores() - 1)

# flog.appender(appender.file(paste0('./log/', format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss"), ".log")))
flog.threshold('trace') %>% invisible

flog.info('Sourcing functions...')
source("src/simulation-study/utils.r")
source("src/simulation-study/generateSamples.r")
source("src/simulation-study/fit.ignorable.r")

source("src/simulation-study/fit.parametric.r")

source("src/simulation-study/fit.hybrid.r")

source("src/simulation-study/fit.class.r")

set.seed(666L)
df1 <- generateSamples(samples = 1, participants = 1000)

# ALL SAMPLES
res <- list()

# res[[1]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% lapply(fit.ignorable) %>% bind_rows

# df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% lapply(fit.parametric) %>% bind_rows

flog.info('Testing hybrid model...')
tic()
res[[1]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% lapply(fit.hybrid) %>% bind_rows
toc()

# flog.info('Fitting models to MAR scenario...')
# res[[1]] <- df1 %>% mutate(y=yMAR, r=rMAR) %>% group_split(sample) %>% mclapply(fit.multiple) %>% bind_rows %>% mutate(scenario='MAR')

# flog.info('Fitting models to MNAR scenario...')
# res[[2]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% mclapply(fit.multiple) %>% bind_rows %>% mutate(scenario='MNAR')

result <- bind_rows(res)


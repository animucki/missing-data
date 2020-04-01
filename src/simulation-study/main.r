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

flog.appender(appender.file(paste0('./log/', format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss"), ".log"))) %>% invisible
flog.threshold('trace') %>% invisible

flog.info('Sourcing functions...')
source("src/simulation-study/utils.r")
source("src/simulation-study/generateSamples.r")
source("src/simulation-study/fit.ignorable.r")

source("src/simulation-study/fit.parametric.r")

source("src/simulation-study/fit.hybrid.r")

# source("src/simulation-study/fit.class.r")

set.seed(666L)
df1 <- generateSamples(samples = 100, participants = 200)

# ALL SAMPLES
res <- list()

# flog.info('Testing class model...')
# tic()
# res[[1]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% lapply(fit.class) %>% bind_rows
# toc()

# flog.info('Fitting models to MAR scenario...')
# res[[1]] <- df1 %>% mutate(y=yMAR, r=rMAR) %>% group_split(sample) %>% mclapply(fit.multiple) %>% bind_rows %>% mutate(scenario='MAR')

flog.info('Fitting models to MNAR scenario...')
res[[2]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% mclapply(fit.multiple) %>% bind_rows %>% mutate(scenario='MNAR')
write.csv2(res[[2]],
           file = './log/outMNAR.csv')

result <- bind_rows(res)


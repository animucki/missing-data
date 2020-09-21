# Run this script to run the simulation.

# clear workspace
rm(list = ls())

library(armspp)
library(futile.logger)
library(lme4)
library(mc2d)
library(numDeriv)
library(parallel)
library(tictoc)
library(tidyverse)

options(mc.cores = parallel::detectCores() - 1)

flog.appender(appender.file(paste0('./log/sim-', format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss"), ".log"))) %>% invisible
flog.threshold('trace') %>% invisible

flog.info('Sourcing functions...')
source("src/common/utils.r")
source("src/simulation-study/generateSamples.r")
source("src/simulation-study/fit.ignorable.r")
source("src/simulation-study/fit.parametric.r")
source("src/simulation-study/fit.hybrid.r")
source("src/simulation-study/fit.class.r")
source("src/simulation-study/fit.npsp.r")
source("src/simulation-study/fit.multiple.r")

set.seed(666L)
df1 <- generateSamples(samples = 100, participants = 200)

# ALL SAMPLES
res <- list()

flog.info('Fitting models to MAR scenario...')
res[[1]] <- df1 %>% mutate(y=yMAR, r=rMAR) %>% group_split(sample) %>% mclapply(fit.multiple, mc.preschedule = F) %>% bind_rows %>% mutate(scenario='MAR')

flog.info('Fitting models to MNAR scenario...')
res[[2]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% mclapply(fit.multiple, mc.preschedule = F) %>% bind_rows %>% mutate(scenario='MNAR')

flog.info('Fitting models to MNAR1 scenario...')
res[[3]] <- df1 %>% mutate(y=yMNAR1, r=rMNAR1) %>% group_split(sample) %>% mclapply(fit.multiple, mc.preschedule = F) %>% bind_rows %>% mutate(scenario='MNAR1')

flog.info('Fitting models to MNAR2 scenario...')
res[[4]] <- df1 %>% mutate(y=yMNAR2, r=rMNAR2) %>% group_split(sample) %>% mclapply(fit.multiple, mc.preschedule = F) %>% bind_rows %>% mutate(scenario='MNAR2')


result <- bind_rows(res)
write.csv2(result, './data/result.csv')

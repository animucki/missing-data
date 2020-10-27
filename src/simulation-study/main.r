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
source("src/simulation-study/fit.tseng.direct.r")
source("src/simulation-study/fit.class.r")
source("src/simulation-study/fit.npsp.r")
source("src/simulation-study/fit.multiple.r")

set.seed(668L)
df1 <- generateSamples(samples = 50, participants = 100)

# ALL SAMPLES
res <- list()

flog.info('Fitting models to MAR scenario...')
tic('mar')
res[[1]] <- df1 %>% mutate(y=yMAR, r=rMAR) %>% group_split(sample) %>% mclapply(fit.multiple, mc.preschedule = F) %>% bind_rows %>% mutate(scenario='MAR')
toc()

flog.info('Fitting models to MNAR scenario...')
tic('mnar')
res[[2]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% mclapply(fit.multiple, mc.preschedule = F) %>% bind_rows %>% mutate(scenario='MNAR')
toc()

flog.info('Fitting models to MNAR3 scenario...')
tic('mnar3')
res[[3]] <- df1 %>% mutate(y=yMNAR3, r=rMNAR3) %>% group_split(sample) %>% mclapply(fit.multiple, mc.preschedule = F) %>% bind_rows %>% mutate(scenario='MNAR3')
toc()

flog.info('Fitting models to MNAR4 scenario...')
tic('mnar4')
res[[4]] <- df1 %>% mutate(y=yMNAR4, r=rMNAR4) %>% group_split(sample) %>% mclapply(fit.multiple, mc.preschedule = F) %>% bind_rows %>% mutate(scenario='MNAR4')
toc()

result <- bind_rows(res)
View(result)
write.csv2(result, './data/result-allmodels-samples-0001-to-0050.csv', row.names = FALSE)

# clear workspace
rm(list = ls())

set.seed(666)

library(dplyr)
library(futile.logger)
library(lme4)
library(fastGHQuad)
library(mvtnorm)
library(numDeriv)
library(stats4)
library(tictoc)
library(tidyr)

flog.info('Sourcing functions...')
source("src/simulation-study/generateSamples.r")
source("src/simulation-study/fit.ignorable.r")
source("src/simulation-study/fit.parametric.r")

flog.appender(appender.file(paste0('./log/', format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss"), ".log")))
flog.threshold('trace')
df <- generateSamples(samples = 1, participants = 10000)

# ALL SAMPLES
res <- list()

flog.info('Fitting ignorable model to MAR scenario...')
# res[[1]] <- df %>% mutate(y=yMAR) %>% group_by(sample) %>% group_modify(fit.ignorable) %>% mutate(model='ignorable', scenario='MAR')

flog.info('Fitting ignorable model to MNAR scenario...')
# res[[2]] <- df %>% mutate(y=yMNAR) %>% group_by(sample) %>% group_modify(fit.ignorable) %>% mutate(model='ignorable', scenario='MNAR')

flog.info('Fitting parametric model to MAR scenario...')
res[[3]] <- df %>% mutate(y=yMAR, r=rMAR) %>% group_by(sample) %>% group_modify(fit.parametric) %>% mutate(model='parametric', scenario='MAR')

flog.info('Fitting parametric model to MNAR scenario...')
# res[[4]] <- df %>% mutate(y=yMNAR, r=rMNAR) %>% group_by(sample) %>% group_modify(fit.parametric) %>% model(model='parametric', scenario='MNAR')

result <- bind_rows(res)


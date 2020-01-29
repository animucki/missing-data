# clear workspace
rm(list = ls())

flog.info('Loading libraries...')
library(dplyr)
library(futile.logger)
library(lme4)
library(mvtnorm)
library(numDeriv)
library(stats4)
library(tidyr)

flog.info('Sourcing functions...')
source("src/simulation-study/generateSamples.r")
source("src/simulation-study/fit.ignorable.r")

flog.threshold('debug')
df <- generateSamples(samples = 10)

# ALL SAMPLES
res <- list()

flog.debug('Fitting ignorable model to MAR scenario...')
res[[1]] <- df %>% mutate(y=yMAR, model='ignorable', scenario='MAR') %>% group_by(sample) %>% group_modify(fit.ignorable)

flog.debug('Fitting ignorable model to MNAR scenario...')
res[[2]] <- df %>% mutate(y=yMNAR, model='ignorable', scenario='MNAR') %>% group_by(sample) %>% group_modify(fit.ignorable)

result <- bind_rows(res)

## ONE SAMPLE
d <- df %>% filter(sample == 1) %>% mutate(y=yMAR)
fit.ignorable(d, 1)
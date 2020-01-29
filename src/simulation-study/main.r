# clear workspace
rm(list = ls())

flog.info('Loading libraries...')
library(dplyr)
library(futile.logger)
library(lme4)
library(fastGHQuad)
library(mvtnorm)
library(numDeriv)
library(stats4)
library(tidyr)

flog.info('Sourcing functions...')
source("src/simulation-study/generateSamples.r")
source("src/simulation-study/fit.ignorable.r")
source("src/simulation-study/fit.parametric.r")

flog.threshold('trace')
df <- generateSamples(samples = 10)

# ALL SAMPLES
res <- list()

flog.info('Fitting ignorable model to MAR scenario...')
# res[[1]] <- df %>% mutate(y=yMAR, model='ignorable', scenario='MAR') %>% group_by(sample) %>% group_modify(fit.ignorable)

flog.info('Fitting ignorable model to MNAR scenario...')
# res[[2]] <- df %>% mutate(y=yMNAR, model='ignorable', scenario='MNAR') %>% group_by(sample) %>% group_modify(fit.ignorable)

flog.info('Fitting parametric model to MAR scenario...')
res[[3]] <- df %>% mutate(y=yMAR, r=rMAR, model='parametric', scenario='MAR') %>% group_by(sample) %>% group_modify(fit.parametric)

flog.info('Fitting parametric model to MNAR scenario...')
# res[[4]] <- df %>% mutate(y=yMNAR, r=rMNAR, model='parametric', scenario='MNAR') %>% group_by(sample) %>% group_modify(fit.parametric)

result <- bind_rows(res)


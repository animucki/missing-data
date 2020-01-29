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

# View(df)

## ONE SAMPLE
d <- df %>% filter(sample == 1) %>% mutate(y=yMAR)
fit.ignorable(d, 1)

# ALL SAMPLES
res <- df %>% mutate(y=yMAR) %>% group_by(sample) %>% group_modify(fit.ignorable)

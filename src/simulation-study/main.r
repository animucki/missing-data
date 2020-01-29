# clear workspace
rm(list = ls())

flog.info('Loading libraries...')
library(dplyr)
library(futile.logger)
library(mvtnorm)
library(numDeriv)
library(stats4)
library(tidyr)

flog.info('Sourcing functions...')
source("src/simulation-study/generateSamples.r")


flog.threshold('info')
df <- generateSamples(samples = 10)

# View(df)

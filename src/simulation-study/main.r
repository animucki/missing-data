# clear workspace
rm(list = ls())

library(dplyr)
library(futile.logger)
library(mvtnorm)
library(numDeriv)
library(stats4)
library(tidyr)

source("src/simulation-study/generateSamples.r")

df <- generateSamples(samples = 10)

View(df)

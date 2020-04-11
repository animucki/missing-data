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

# options(mc.cores = parallel::detectCores() - 1)

# flog.appender(appender.file(paste0('./log/emp-', format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss"), ".log"))) %>% invisible
flog.threshold('trace') %>% invisible

flog.info('Sourcing functions...')
source("src/common/utils.r")


flog.info('Loading data...')
data <- read.csv('./data/scleroderma.csv')

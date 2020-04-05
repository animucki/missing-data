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

source("src/simulation-study/fit.class.r")

set.seed(666L)
df1 <- generateSamples(samples = 200, participants = 200)

# ALL SAMPLES
res <- list()

# flog.info('Testing class model...')
# tic()
# res[[1]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% lapply(fit.class) %>% bind_rows
# toc()

flog.info('Fitting models to MAR scenario...')
res[[1]] <- df1 %>% mutate(y=yMAR, r=rMAR) %>% group_split(sample) %>% mclapply(fit.multiple) %>% bind_rows %>% mutate(scenario='MAR')
# write.csv2(res[[1]], file = '/home/bartosz/Dropbox/outSurvMAR.csv')

flog.info('Fitting models to MNAR scenario...')
res[[2]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% mclapply(fit.multiple) %>% bind_rows %>% mutate(scenario='MNAR')
# write.csv2(res[[2]], file = '/home/bartosz/Dropbox/outSurvMNAR.csv')

# result <- bind_rows(res)

stop('It\'s all over now')

result <- read.csv2('./result.csv', stringsAsFactors = F, row.names = 1)
result <- result %>%
  mutate(bias_intercept = -1.2 - intercept,
         bias_time = 0.5 - time,
         bias_treatment = -1.5 - treatment,
         bias_sigma.b = 0.5 - sigma.b,
         bias_sigma = sqrt(0.5) - sigma,
         coverage_intercept = abs(bias_intercept) <= se.intercept * qnorm(0.975),
         coverage_time = abs(bias_time) <= se.time * qnorm(0.975),
         coverage_treatment = abs(bias_treatment) <= se.treatment * qnorm(0.975),
         coverage_sigma.b = abs(bias_sigma.b) <= se.sigma.b * qnorm(0.975),
         coverage_sigma = abs(bias_sigma) <= se.sigma * qnorm(0.975),
         length_intercept = 2 * se.intercept * qnorm(0.975),
         length_time = 2 * se.time * qnorm(0.975),
         length_treatment = 2 * se.treatment * qnorm(0.975),
         length_sigma.b = 2 * se.sigma.b * qnorm(0.975),
         length_sigma = 2 * se.sigma * qnorm(0.975)) %>%
  group_by(scenario, model) %>%
  summarize(avg_bias_intercept = mean(bias_intercept),
            avg_bias_time = mean(bias_time),
            avg_bias_treatment = mean(bias_treatment),
            avg_bias_sigma.b = mean(bias_sigma.b),
            avg_bias_sigma = mean(bias_sigma),
            mcse_bias_intercept = sd(bias_intercept),
            mcse_bias_time = sd(bias_time),
            mcse_bias_treatment = sd(bias_treatment),
            mcse_bias_sigma.b = sd(bias_sigma.b),
            mcse_bias_sigma = sd(bias_sigma),
            coverage_intercept = mean(coverage_intercept),
            coverage_time = mean(coverage_time),
            coverage_treatment = mean(coverage_treatment),
            coverage_sigma.b = mean(coverage_sigma.b),
            coverage_sigma = mean(coverage_sigma),
            length_intercept = mean(length_intercept),
            length_time = mean(length_time),
            length_treatment = mean(length_treatment),
            length_sigma.b = mean(length_sigma.b),
            length_sigma = mean(length_sigma))

# View(result)

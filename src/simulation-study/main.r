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

#this is SPECIAL CODE!
resPartial <- read.csv2('./data/result.csv', stringsAsFactors = F, row.names = 1) %>%
  filter(model == "class")

resJobs <- expand.grid(scenario = c("MAR","MNAR"), sample = 1:100, stringsAsFactors = F)
resJobs <- anti_join(resJobs, resPartial)
jobs <- split(resJobs, seq(nrow(resJobs)))

res <- jobs %>% mclapply(function (row) {
  if(row$scenario == "MAR") {
    df1 %>% mutate(y=yMAR, r=rMAR) %>% filter(sample == row$sample) %>% fit.class %>% mutate(scenario='MAR')
  } else {
    df1 %>% mutate(y=yMNAR, r=rMNAR) %>% filter(sample == row$sample) %>% fit.class %>% mutate(scenario='MNAR')
  }}, mc.preschedule = F)

# res[[1]] <- df1 %>% mutate(y=yMAR, r=rMAR) %>% group_split(sample) %>% mclapply(fit.multiple, mc.preschedule = F) %>% bind_rows %>% mutate(scenario='MAR')
# res[[2]] <- df1 %>% mutate(y=yMNAR, r=rMNAR) %>% group_split(sample) %>% mclapply(fit.multiple, mc.preschedule = F) %>% bind_rows %>% mutate(scenario='MNAR')

result <- bind_rows(res)
write.csv2(result, './data/result_p2.csv')

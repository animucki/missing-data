rm(list=ls())
#library(dplyr)
#library(ggplot2)
library(gridExtra)
library(mc2d)
library(tidyverse)

n <- 30

set.seed(666L)
nioMAR <- 0
nioMCAR <- 1
while (nioMAR != nioMCAR) {
  y <- rmultinormal(n, c(125,125), 25^2*c(1, 0.6, 0.6, 1))
  nioMAR <- sum(y[,1] > 140)
  nioMCAR <- sum(y[,2] > 140)
}

d <- data.frame(y1 = y[,1],
                y2 = y[,2],
                y2MCAR = NA_real_) %>%
  mutate(y2MAR = case_when(
    y1 > 140 ~ y2,
    TRUE ~ NA_real_
  ),
         y2MNAR = case_when(
    y2 > 140 ~ y2,
    TRUE ~ NA_real_
  ))

indices_obs_MCAR <- sample.int(n,nioMAR)
d$y2MCAR[indices_obs_MCAR] <- d$y2[indices_obs_MCAR]

subplots <- list()
subplots[[1]] <- ggplot(d, aes(x=y1, y=y2)) +
  geom_point(shape=1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method=lm, se=FALSE) +
  xlim(range(d$y1)) +
  ylim(range(d$y2)) +
  labs(title="(a) Complete data",
       x=expression(Y[1]), y=expression(Y[2]))
subplots[[2]] <- ggplot(d, aes(x=y1, y=y2MCAR)) +
  geom_point(shape=1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method=lm, se=FALSE) +
  xlim(range(d$y1)) +
  ylim(range(d$y2)) +
  labs(title="(b) Randomly sampled points",
       x=expression(Y[1]), y=expression(Y[2]))
subplots[[3]] <- ggplot(d, aes(x=y1, y=y2MAR)) +
  geom_point(shape=1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method=lm, se=FALSE) +
  xlim(range(d$y1)) +
  ylim(range(d$y2)) +
  labs(title=expression("(c) Missing if"~Y[1]<140),
       x=expression(Y[1]), y=expression(Y[2]))
subplots[[4]] <- ggplot(d, aes(x=y1, y=y2MNAR)) +
  geom_point(shape=1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method=lm, se=FALSE) +
  xlim(range(d$y1)) +
  ylim(range(d$y2)) +
  labs(title=expression("(d) Missing if"~Y[2]<140),
       x=expression(Y[1]), y=expression(Y[2]))
(plotArranged <- grid.arrange(grobs = subplots, top = "Bivariate example of missing data mechanisms"))
ggsave("./plots/introductionPlot.pdf", plot = plotArranged,
       width = 10,
       height = 7.5,
       units = "in",
       device = cairo_pdf)


dRes <- NULL
for (yCol in 1:(ncol(d)-1)) {
  mu_y <- mean(d[,1+yCol],na.rm = T)
  se_y <- sd(d[,1+yCol],na.rm = T)/sqrt(sum(!is.na(d[,1+yCol])))

  m <- lm(d[,1+yCol] ~ d$y1)

  print(rbind(mu_y + c(0, -se_y, +se_y)*qt(0.975,n-1), cbind(coef(m),confint(m))))

  dRes <- rbind(dRes, c(mu_y, coef(m), se_y, summary(m)$coefficients[,2], use.names = F))
}
colnames(dRes) <- c("est.mu2", "est.beta1", "est.beta2", "se.mu2", "se.beta1", "se.beta2")
dRes <- as.data.frame(dRes)
dRes$scenario <- c("(a)","(b)","(c)","(d)")
dRes1 <- dRes %>%
  select(-starts_with("se")) %>%
  pivot_longer(starts_with("est"), values_to = "est") %>%
  mutate(estimand = substring(name, 5)) %>%
  select(-name)
dRes2 <- dRes %>%
  select(-starts_with("est")) %>%
  pivot_longer(starts_with("se"), values_to = "se") %>%
  mutate(estimand = substring(name, 4)) %>%
  select(-name)
dRes <- full_join(dRes1, dRes2)

pRes <- list()
pRes[[1]] <- ggplot(dRes %>% filter(estimand=='mu2'), aes(x=scenario, y=est)) +
  geom_pointrange(aes(ymin=est-se*qt(0.975,n-1), ymax=est+se*qt(0.975,n-1))) +
  geom_hline(yintercept = 125, linetype="dotted") +
  ylab(expression(widehat(mu[2]))) +
  ggtitle(expression("Marginal mean of"~Y[2]))

pRes[[2]] <- ggplot(dRes %>% filter(estimand=='beta1'), aes(x=scenario, y=est)) +
  geom_pointrange(aes(ymin=est-se*qt(0.975,n-1), ymax=est+se*qt(0.975,n-1))) +
  geom_hline(yintercept = 0, linetype="dotted") +
  ylab(expression(widehat(beta[1]))) +
  ggtitle("Regression intercept")

pRes[[3]] <- ggplot(dRes %>% filter(estimand=='beta2'), aes(x=scenario, y=est)) +
  geom_pointrange(aes(ymin=est-se*qt(0.975,n-1), ymax=est+se*qt(0.975,n-1))) +
  geom_hline(yintercept = 1, linetype="dotted") +
  ylab(expression(widehat(beta[2]))) +
  ggtitle("Regression slope")

(pResTotal <- grid.arrange(grobs = pRes, widths=c(1,1,1), top = "Estimated parameters in blood pressure data"))
ggsave("./plots/introductionPlotEstimates.pdf", plot = pResTotal,
       width = 10,
       height = 3.5,
       units = "in",
       device = cairo_pdf)
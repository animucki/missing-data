library(ggplot2)
library(grid)
library(gridExtra)
library(tidyverse)

# Summary of measures for table
result <- read.csv2('./data/result.csv', stringsAsFactors = F, row.names = 1)
resultSummary <- result %>%
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
            empse_intercept = sd(bias_intercept),
            empse_time = sd(bias_time),
            empse_treatment = sd(bias_treatment),
            empse_sigma.b = sd(bias_sigma.b),
            empse_sigma = sd(bias_sigma),
            coverage_intercept = 100 * mean(coverage_intercept),
            coverage_time = 100 * mean(coverage_time),
            coverage_treatment = 100 * mean(coverage_treatment),
            coverage_sigma.b = 100 * mean(coverage_sigma.b),
            coverage_sigma = 100 * mean(coverage_sigma),
            length_intercept = mean(length_intercept),
            length_time = mean(length_time),
            length_treatment = mean(length_treatment),
            length_sigma.b = mean(length_sigma.b),
            length_sigma = mean(length_sigma))

resultTable <- list()

resultTable[[1]] <- t(resultSummary[,-c(1,2)])[1:5,]
colnames(resultTable[[1]]) <- resultSummary %>% mutate(tmp = paste(scenario,model)) %>% pull(tmp)
resultTable[[1]] <- format(resultTable[[1]], digits = 1, scientific = F, trim = T)

resultTable[[2]] <- t(resultSummary[,-c(1,2)])[6:10,]
colnames(resultTable[[2]]) <- resultSummary %>% mutate(tmp = paste(scenario,model)) %>% pull(tmp)
resultTable[[2]] <- format(resultTable[[2]], digits = 2, scientific = F, trim = T)
resultTable[[2]] <- matrix(paste0("(",resultTable[[2]],")"),5,10)

resultTable[[3]] <- t(resultSummary[,-c(1,2)])[11:15,]
colnames(resultTable[[3]]) <- resultSummary %>% mutate(tmp = paste(scenario,model)) %>% pull(tmp)
resultTable[[3]] <- format(resultTable[[3]], nsmall = 1, trim = T)

resultTable[[4]] <- t(resultSummary[,-c(1,2)])[16:20,]
colnames(resultTable[[4]]) <- resultSummary %>% mutate(tmp = paste(scenario,model)) %>% pull(tmp)
resultTable[[4]] <- format(resultTable[[4]], digits = 1, trim = T)
resultTable[[4]] <- matrix(paste0("(",resultTable[[4]],")"),5,10)

order <- as.vector(matrix(1:10, nrow = 2, byrow = T))
order <- c(order, 10+order)
resultTableCombined <- do.call(rbind, resultTable)[order,]
colOrder <- c(3,4,2,1,5)
resultTableCombined[,colOrder]
resultTableCombined[,5+colOrder]

problematicCells <- abs(t(resultSummary[,-c(1,2)])[1:5,]) > t(resultSummary[,-c(1,2)])[6:10,]
problematicCells[,colOrder]
problematicCells[,5+colOrder]

#Type 1 plots: univariate plots to assess the overall distribution and outliers

result2 <- result
result2$scenario <- factor(result2$scenario, levels = c("MAR","MNAR"), ordered = T)
result2$model <- factor(result2$model, levels = c("ignorable","spm","hybrid","class","spsp"), ordered = T)

estimands <- c('intercept', 'time', 'treatment', 'sigma.b', 'sigma')
estimandsNice <- c(expression(widehat(beta[1])),
                   expression(widehat(beta[2])),
                   expression(widehat(beta[3])),
                   expression(widehat(sigma[b])),
                   expression(widehat(sigma)))
seNice <- c(expression(widehat(SE) * (widehat(beta[1]))),
            expression(widehat(SE) * (widehat(beta[2]))),
            expression(widehat(SE) * (widehat(beta[3]))),
            expression(widehat(SE) * (widehat(sigma[b]))),
            expression(widehat(SE) * (widehat(sigma))))
trueValues <- c(-1.2, 0.5, -1.5, 0.5, sqrt(0.5))

plots1_estimand <- list()
for(i in seq_along(estimands)) {
  plots1_estimand[[i]] <-
    ggplot(result2, aes_string(x="model", y=estimands[i], fill="scenario")) +
      geom_boxplot(outlier.shape = 1) +
      geom_hline(yintercept = trueValues[i], linetype="dotted") +
      scale_fill_brewer(palette = "Paired") +
      ylab(estimandsNice[i]) +
      # ggtitle(paste("Distribution of estimated",estimands[i])) +
      theme(axis.title.x = element_blank())

  if(i<5) plots1_estimand[[i]] <- plots1_estimand[[i]] + theme(legend.title = element_blank(), legend.position = "none")
}
# plots1_estimand

plots1_se <- list()
for(i in seq_along(estimands)) {
  plots1_se[[i]] <-
    ggplot(result2, aes_string(x="model", y=paste0("se.",estimands[i]), fill="scenario")) +
      geom_boxplot(outlier.shape = 1) +
      scale_fill_brewer(palette = "Paired") +
      ylab(seNice[i]) +
      # ggtitle("Distribution of estimated SE of",estimands[i]) +
      theme(axis.title.x = element_blank())

  if(i<5) plots1_se[[i]] <- plots1_se[[i]] + theme(legend.title = element_blank(), legend.position = "none")
}
# plots1_se

#Export all type 1 plots
for(i in seq_along(estimands)) {
  if(i==5) {
    ggsave(paste0("./plots/plot1_est_", estimands[i],".pdf"),
           plot = plots1_estimand[[i]],
           device = cairo_pdf,
           width = 6,
           height = 3,
           units = "in")
    ggsave(paste0("./plots/plot1_se_", estimands[i],".pdf"),
           plot = plots1_se[[i]],
           device = cairo_pdf,
           width = 6,
           height = 3,
           units = "in")
  } else {
    ggsave(paste0("./plots/plot1_est_", estimands[i],".pdf"),
           plot = plots1_estimand[[i]],
           device = cairo_pdf,
           width = 5,
           height = 3,
           units = "in")
    ggsave(paste0("./plots/plot1_se_", estimands[i],".pdf"),
           plot = plots1_se[[i]],
           device = cairo_pdf,
           width = 5,
           height = 3,
           units = "in")
  }
}

#Type 2 plots: bivariate plots of est.se vs est. for each model/scenario combination
combinations <- expand.grid(estimand=estimands, scenario=levels(result2$scenario),model=levels(result2$model))[,3:1] %>%
  filter(!(model=='spsp' & estimand =='sigma.b'))

plots2 <- list()
for(i in seq_len(nrow(combinations))) {
  current_estNice <- estimandsNice[which(estimands==combinations[i,"estimand"])]
  current_seNice <- seNice[which(estimands==combinations[i,"estimand"])]

  plots2[[i]] <- ggplot(result2 %>% filter(model==as.character(combinations[i,"model"]), scenario==as.character(combinations[i,"scenario"])),
                        aes_string(x=as.character(combinations[i,"estimand"]), y=paste0("se.",combinations[i,"estimand"]))) +
    geom_density2d() +
    geom_point() +
    xlab(current_estNice) +
    ylab(current_seNice) #+
    # ggtitle( paste0("Est. SE vs estimand for ", combinations[i,"model"], " model under ", combinations[i,"scenario"]) )
}
# plots2

#Export the type 2 plots
for(i in seq_len(nrow(combinations))) {
  ggsave(paste("./plots/plot2", combinations[i,1], combinations[i,2], combinations[i,3], ".pdf", sep="_"),
         width = 2*2.5,
         height = 2*1.875,
         units = "in",
         plot = plots2[[i]],
         device = cairo_pdf
  )
}

#Type 3 plots: bivariate plots of estimand/SE between every pair of models
plots3 <- list()
combinations3 <- expand.grid(estimand=estimands, scenario=levels(result2$scenario))
for(i in seq_len(nrow(combinations3))) {
  currentData <- result2 %>% filter(as.character(scenario)==combinations3[i,"scenario"]) %>%
    select("sample","model", as.character(combinations3[i,"estimand"]), paste0("se.",combinations3[i,"estimand"]))
  models <- levels(currentData$model)

  stopifnot(length(models) > 1)

  pairs <- expand.grid(modelX=models,modelY=models) %>%
    filter(as.integer(modelX) > as.integer(modelY))

  subplots_est <- list()
  for (j in seq_len(nrow(pairs))) {
    subData <- currentData %>% filter(model %in% unlist(pairs[j,])) %>%
      pivot_wider(id_cols = sample, names_from = model, values_from = 3)

    subplots_est[[j]] <- ggplot(subData, aes_string(x=as.character(pairs[j,1]),
                                                    y=as.character(pairs[j,2]) )) +
      geom_point() +
      geom_abline(slope=1, linetype="dotted") +
      theme(axis.title = element_blank(),
            axis.text.x = element_text(angle = 45))
  }

  subplots_se <- list()
  for (j in seq_len(nrow(pairs))) {
    subData <- currentData %>% filter(model %in% unlist(pairs[j,])) %>%
      pivot_wider(id_cols = sample, names_from = model, values_from = 4)

    subplots_se[[j]] <- ggplot(subData, aes_string(x=as.character(pairs[j,2]),
                                                   y=as.character(pairs[j,1]) )) +
      geom_point() +
      geom_abline(slope=1, linetype="dotted") +
      theme(axis.title = element_blank(),
            axis.text.x = element_text(angle = 45))
  }

  modelNames <- lapply(models, function(t) textGrob(t))
  subplots <- c(subplots_est, subplots_se, modelNames)

  plots3[[i]] <- grid.arrange(grobs = subplots,
                              layout_matrix = rbind(c(21, 1, 2, 3, 4),
                                                    c(11,22, 5, 6, 7),
                                                    c(12,15,23, 8, 9),
                                                    c(13,16,18,24,10),
                                                    c(14,17,19,20,25)),
                              top=textGrob(label = paste0("Between-model plots for ", as.character(combinations3[i,"estimand"]),
                                                          " under ", as.character(combinations3[i,"scenario"])) ))
}

#export the type 3 plots
for (i in seq_len(nrow(combinations3))) {
  ggsave(paste("./plots/plot3", combinations3[i,1], combinations3[i,2], ".pdf", sep="_"),
         width = 10,
         height = 10,
         units = "in",
         plot = plots3[[i]],
         device = cairo_pdf
  )
}

#' Calculate the multinomial logit prediction function
softmax <- function(x) {
  exp(x) / sum(exp(x))
}

#' Fit selected models for one sample
fit.multiple <- function(d, models = c('ignorable','spm','spm+class','class')) {
  res <- list()
  if('ignorable' %in% models) res[[1]] <- d %>% fit.ignorable  %>% mutate(model='ignorable')
  if('spm'       %in% models) res[[2]] <- d %>% fit.parametric %>% mutate(model='spm')
  if('spm+class' %in% models) res[[3]] <- d %>% fit.hybrid     %>% mutate(model='spm+class')
  # if('class'     %in% models) res[[4]] <- d %>% fit.class      %>% mutate(model='class')
  return(bind_rows(res))
}
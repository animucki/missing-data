#' Fit selected models for one sample
fit.multiple <- function(d, models = c('ignorable','spm','tseng','class','spsp')) {
  res <- list()
  #if('ignorable' %in% models) res[[1]] <- d %>% fit.ignorable  %>% mutate(model='ignorable')
  #if('heckman'   %in% models) res[[2]] <- d %>% fit.heckman    %>% mutate(model='heckman')
  #if('spm'       %in% models) res[[3]] <- d %>% fit.parametric %>% mutate(model='spm')
  if('tseng'     %in% models) res[[4]] <- d %>% fit.tseng      %>% mutate(model='tseng')
  #if('class'     %in% models) res[[5]] <- d %>% fit.class      %>% mutate(model='class')
  #if('spsp'      %in% models) res[[6]] <- d %>% fit.npsp       %>% mutate(model='spsp')

  res <- bind_rows(res)

  return(res)
}
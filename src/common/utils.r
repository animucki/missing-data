#' Calculate the multinomial logit prediction function
softmax <- function(x) {
  t <- exp(x) / sum(exp(x))
  if(any(!is.finite(t))) {
    return(as.numeric(x==max(x)))
  } else {
    return(t)
  }
}

#' matrix trace
tr <- function(x) {
  sum(diag(x))
}
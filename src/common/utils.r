#' Calculate the multinomial logit prediction function
softmax <- function(x) {
  exp(x) / sum(exp(x))
}

#' matrix trace
tr <- function(x) {
  sum(diag(x))
}
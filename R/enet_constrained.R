#' @keywords internal
expit <- function(x) exp(x) / (1 + exp(x))

#' @keywords internal
logit <- function(x) log(x) - log(1 - x)

#' @keywords internal
norm_l2 <- function(beta) sqrt(sum(beta ^ 2))

#' @keywords internal
norm_l1 <- function(beta) sum(abs(beta))


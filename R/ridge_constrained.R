#' @keywords internal
expit <- function(x) exp(x) / (1 + exp(x))
cross_entropy <- function(beta, y, x) {
  y_is_one <- y == 1
  l_one <- -log(expit(x %*% beta))
  l_zero <- -log(1 - expit(x %*% beta))
  l <- l_zero
  l[y_is_one] <- l_one[y_is_one]
  return(mean(l))
}
norm_l2 <- function(beta) sqrt(sum(beta ^ 2))


#' logistic ridge regression (constrained form)
#'
#' @param Y vector of binary outcome
#' @param X matrix of predictor
#' @param beta_init initial value for slope
#' @param l2_norm_max upper constraint on l2 norm of slope
#'
#' @return the fitted slope
#' @export
fit_ridge_constrained <- function(Y, X, beta_init, l2_norm_max) {
  ridge_fit <- Rsolnp::solnp(
    pars = beta_init,
    fun = function(b) cross_entropy(b, y = Y, x = X),
    ineqfun = norm_l2,
    ineqLB = 0,
    ineqUB = l2_norm_max,
    control = list(trace = 0)
  )
  sol <- ridge_fit$pars
  return(sol)
}

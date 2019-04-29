#' @keywords internal
expit <- function(x) exp(x) / (1 + exp(x))

#' @keywords internal
logit <- function(x) log(x) - log(1 - x)

#' @keywords internal
norm_l2 <- function(beta) sqrt(sum(beta ^ 2))

#' @keywords internal
norm_l1 <- function(beta) sum(abs(beta))

# #' @keywords internal
# cross_entropy <- function(beta, y, x, offset) {
#   y_is_one <- y == 1
#   l_one <- -log(expit(x %*% beta))
#   l_zero <- -log(1 - expit(x %*% beta + offset))
#   l <- l_zero
#   l[y_is_one] <- l_one[y_is_one]
#   return(mean(l))
# }


# #' logistic elastic net regression (constrained form)
# #'
# #' @param Y vector of binary outcome
# #' @param X matrix of predictor
# #' @param beta_init initial value for slope
# #' @param norm_max upper constraint on l2/l1 norm of slope
# #' @param offset a vector of constant offset in the linear term
# #' @param type fit ridge regression ("l2") or lasso ("l1")
# #'
# #' @return the fitted slope
# #' @importFrom Rsolnp solnp
# #' @export
# fit_enet_constrained <- function(Y, X, beta_init, norm_max, offset = NULL, type = "l2") {
#   if (is.null(offset)) offset <- rep(0, length(Y))
#   if (type == "l2") constraint_func <- norm_l2
#   if (type == "l1") constraint_func <- norm_l1
#   ridge_fit <- Rsolnp::solnp(
#     pars = beta_init,
#     fun = function(b) cross_entropy(b, y = Y, x = X, offset = offset),
#     ineqfun = constraint_func,
#     ineqLB = 0,
#     ineqUB = norm_max,
#     control = list(trace = 0)
#   )
#   sol <- ridge_fit$pars
#   return(sol)
# }

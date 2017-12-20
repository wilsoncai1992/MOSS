#' Plot the survival curve estimator
#'
#' @param onestepfit object returned by surv_onestep or surv_onestep_complete
#' @param col line color
#' @param add whether to add to existing plot
#' @param ...
#'
#' @return NA
#' @export
#'
#' @examples
# plot.surv_onestep <- function(onestepfit, col = 'green', add = FALSE, ...) {
#   step_curve <- stepfun(x = onestepfit$T.uniq, y = c(1, onestepfit$Psi.hat))
#   curve(step_curve, from = 0, to = max(onestepfit$T.uniq), add = add, col = col, ...)
# }

#' Plot the survival curve estimator
#'
#' @param obj
#' @param add
#' @param ...
#'
#' @return NA
#' @export
#'
#' @examples
plot.surv_survtmle <- function(obj, add = FALSE, ...) {
  step_curve <- stepfun(x = obj$T.uniq, y = c(1, obj$s_vec), right = TRUE)
  curve(step_curve, from = 0, to = max(obj$T.uniq), add = add, ...)
}


#' Title
#'
#' @param fit_obj
#' @param q
#' @param add
#' @param col
#' @param ...
#'
#' @return NA
#' @export
#'
#' @examples
plot_CI <- function(fit_obj, q = 0.95, add = FALSE, col = 'black', ...) {
  sd_CI <- sqrt(fit_obj$var)
  upper <- fit_obj$Psi.hat + qnorm(p = (1-q)/2, lower.tail = FALSE) * sd_CI
  lower <- fit_obj$Psi.hat - qnorm(p = (1-q)/2, lower.tail = FALSE) * sd_CI

  # ad-hoc thresholding between (0,1)
  upper[upper > 1] <- 1
  upper[upper < 0] <- 0
  lower[lower > 1] <- 1
  lower[lower < 0] <- 0

  step_curve_upper <- stepfun(x = fit_obj$T.uniq, y = c(1, upper))
  curve(step_curve_upper, from = 0, to = max(fit_obj$T.uniq), add = add, col = col, ...)

  step_curve_lower <- stepfun(x = fit_obj$T.uniq, y = c(1, lower))
  curve(step_curve_lower, from = 0, to = max(fit_obj$T.uniq), add = add, col = col, ...)

  return(list(upper = upper, lower = lower))
}


#' Title
#'
#' @param onestepfit
#' @param add
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_initial_G <- function(onestepfit, add = FALSE, ...) {
  G_hat <- colMeans(onestepfit$initial_fit$G.hat.t$out_censor_full)
  T.uniq <- onestepfit$T.uniq
  if(add){
    lines(y = G_hat, x = 1:max(T.uniq), ...)
  }else{
    plot(y = G_hat, x = 1:max(T.uniq), type = 'l', ...)
  }
}


#' Title
#'
#' @param onestepfit
#' @param add
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_initial_S <- function(onestepfit, add = FALSE, ...) {
  S_hat <- colMeans(onestepfit$initial_fit$Qn.A1.t)
  T.uniq <- onestepfit$T.uniq
  if(add){
    lines(S_hat ~ T.uniq, ...)
  }else{
    plot(S_hat ~ T.uniq, type = 'l', ...)
  }
}


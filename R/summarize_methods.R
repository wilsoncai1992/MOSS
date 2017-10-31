#' summarize propensity score weights from surv_onestep object
#'
#'
#' @param onestep_curve
#'
#' @return
#' @export
#'
#' @examples
#' # NA
summarize_weights <- function(onestep_curve) {
  out <- list(
    summary(as.vector(as.matrix(onestep_curve$initial_fit$G.hat.t$out_censor))),
    summary(onestep_curve$initial_fit$g.fitted),
    summary(as.vector(as.matrix(onestep_curve$initial_fit$G.hat.t$out_censor * onestep_curve$initial_fit$g.fitted))),
    summary(1/as.vector(as.matrix(onestep_curve$initial_fit$G.hat.t$out_censor * onestep_curve$initial_fit$g.fitted)))
  )
  return(out)
}

require("R6")

#' run ipcw based on initial fit
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods
#' @format \code{\link{R6Class}} object.
#' @examples
#' \donttest{
#'    ipcw$new(A, T_tilde, Delta, density_failure, density_censor, g1W, A_intervene)
#' }
#' @field A vector of treatment
#' @field T_tilde vector of last follow up time
#' @field Delta vector of censoring indicator
#' @field density_failure survival_curve object of predicted counterfactual
#'  survival curve
#' @field density_censor survival_curve object of predicted counterfactual
#'  failure event survival curve
#' @field g1W propensity score
#' @field A_intervene the intervention of interest
#' @field k_grid vector of interested time points
#' @section Methods:
#' fit fit one time point
#' @export
ipcw <- R6Class("ipcw",
  public = list(
    A = NULL,
    T_tilde = NULL,
    Delta = NULL,
    density_failure = NULL,
    density_censor = NULL,
    g1W = NULL,
    A_intervene = NULL,
    initialize = function(
      A, T_tilde, Delta, density_failure, density_censor, g1W, A_intervene
    ) {
      self$A <- A
      self$T_tilde <- T_tilde
      self$Delta <- Delta
      self$density_failure <- density_failure
      self$density_censor <- density_censor
      self$g1W <- g1W
      self$A_intervene <- A_intervene
      return(self)
    },
    fit = function(k) {
      v1 <- self$density_failure$survival[, k]
      v2 <- self$Delta / self$density_censor$survival[
        cbind(1:self$density_censor$n(), self$T_tilde)
      ]
      if (self$A_intervene == 1) g <- self$g1W  else g <- 1 - self$g1W
      v3 <- as.numeric(self$A == self$A_intervene) / g
      return(mean(v1 * v2 * v3))
    }
  )
)

#' run ee based on initial fit
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods
#' @format \code{\link{R6Class}} object.
#' @examples
#' \donttest{
#'    ee$new(A, T_tilde, Delta, density_failure, density_censor, g1W, A_intervene)
#' }
#' @field A vector of treatment
#' @field T_tilde vector of last follow up time
#' @field Delta vector of censoring indicator
#' @field density_failure survival_curve object of predicted counterfactual
#'  survival curve
#' @field density_censor survival_curve object of predicted counterfactual
#'  failure event survival curve
#' @field g1W propensity score
#' @field A_intervene the intervention of interest
#' @field k_grid vector of interested time points
#' @section Methods:
#' fit fit one time point
#' @export
ee <- R6Class("ee",
  public = list(
    A = NULL,
    T_tilde = NULL,
    Delta = NULL,
    density_failure = NULL,
    density_censor = NULL,
    g1W = NULL,
    A_intervene = NULL,
    ipcw = NULL,
    initialize = function(
      A, T_tilde, Delta, density_failure, density_censor, g1W, A_intervene
    ) {
      self$A <- A
      self$T_tilde <- T_tilde
      self$Delta <- Delta
      self$density_failure <- density_failure
      self$density_censor <- density_censor
      self$g1W <- g1W
      self$A_intervene <- A_intervene
      self$ipcw <- ipcw$new(
        A = A,
        T_tilde = T_tilde,
        Delta = Delta,
        density_failure = density_failure,
        density_censor = density_censor,
        g1W = g1W,
        A_intervene = A_intervene
      )
      return(self)
    },
    fit = function(k) {
      psi_n <- self$ipcw$fit(k = k)
      eic_fit <- eic$new(
        A = self$A,
        T_tilde = self$T_tilde,
        Delta = self$Delta,
        density_failure = self$density_failure,
        density_censor = self$density_censor,
        g1W = self$g1W,
        # eic function expect an entire survival curve psi. we only have psi for
        # a single time point. hack around
        psi = rep(psi_n, max(self$T_tilde)),
        A_intervene = self$A_intervene
      )$one_t(k = k)
      return(psi_n + mean(eic_fit))
    }
  )
)

#' repeat ipcw/ee for all time points
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods
#' @format \code{\link{R6Class}} object.
#' @examples
#' \donttest{
#'    repeat_t_grid$new(method = ipcw,
#'      A,
#'      T_tilde,
#'      Delta,
#'      density_failure,
#'      density_censor,
#'      g1W,
#'      A_intervene
#'    )
#'    repeat_t_grid$new(method = ee,
#'      A,
#'      T_tilde,
#'      Delta,
#'      density_failure,
#'      density_censor,
#'      g1W,
#'      A_intervene
#'    )
#' }
#' @field method either ipcw or ee class
#' @field A vector of treatment
#' @field T_tilde vector of last follow up time
#' @field Delta vector of censoring indicator
#' @field density_failure survival_curve object of predicted counterfactual
#'  survival curve
#' @field density_censor survival_curve object of predicted counterfactual
#'  failure event survival curve
#' @field g1W propensity score
#' @field A_intervene the intervention of interest
#' @field k_grid vector of interested time points
#' @section Methods:
#' fit fit all time points in `k_grid`
#' @export
repeat_t_grid <- R6Class("repeat_t_grid",
  public = list(
    method_fit = NULL,
    initialize = function(
      method, A, T_tilde, Delta, density_failure, density_censor, g1W, A_intervene
    ) {
      self$method_fit <- method$new(
        A = A,
        T_tilde = T_tilde,
        Delta = Delta,
        density_failure = density_failure,
        density_censor = density_censor,
        g1W = g1W,
        A_intervene = A_intervene
      )
      return(self)
    },
    fit = function(k_grid) {
      psi_n <- c()
      for (k in k_grid) {
        psi_n <- c(psi_n, self$method_fit$fit(k = k))
      }
      return(psi_n)
    }
  )
)

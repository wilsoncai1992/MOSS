require("R6")

#' compute eic under the likelihood fit for each observation
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods
#' @format \code{\link{R6Class}} object.
#' @field A vector of treatment
#' @field T_tilde vector of last follow up time
#' @field Delta vector of censoring indicator
#' @field density_failure survival_curve object of predicted counterfactual
#'  survival curve
#' @field density_censor survival_curve object of predicted counterfactual
#'  failure event survival curve
#' @field g1W propensity score
#' @field psi a vector of target parameter estimate
#' @field A_intervene the intervention of interest
#' @field k_grid vector of interested time points
#' @section Methods:
#' one_t compute a vector of EIC
#' all_t compute a matrix of EIC for all time points
#' clever_covariate compute the clever covariate for one time point
#' @export
eic <- R6Class("eic",
  public = list(
    A = NULL,
    T_tilde = NULL,
    Delta = NULL,
    density_failure = NULL,
    density_censor = NULL,
    g1W = NULL,
    psi = NULL,
    A_intervene = NULL,
    initialize = function(
      A, T_tilde, Delta, density_failure, density_censor, g1W, psi, A_intervene
    ) {
      self$A <- A
      self$T_tilde <- T_tilde
      self$Delta <- Delta
      self$density_failure <- density_failure
      self$density_censor <- density_censor
      self$g1W <- g1W
      self$psi <- psi
      self$A_intervene <- A_intervene
      return(self)
    },
    one_t = function(k) {
      if (self$A_intervene == 1) g <- self$g1W  else g <- 1 - self$g1W
      part1_sum <- rep(0, length(g))
      for (t in 1:k) {
        h <- -as.numeric(self$A == self$A_intervene) / g /
          self$density_censor$survival[, t] *
          self$density_failure$survival[, k] / self$density_failure$survival[, t]
        part1 <-  h * (
          as.numeric(self$T_tilde == t & self$Delta == 1) -
          as.numeric(self$T_tilde >= t) * self$density_failure$hazard[, t]
        )
        part1_sum <- part1_sum + part1
      }
      part2 <- self$density_failure$survival[, k] - self$psi[k]
      return(part1_sum + part2)
    },
    all_t = function(k_grid) {
      # naive way to compute for all t
      eic_all <- list()
      for (k in k_grid) {
        eic_all <- c(eic_all, list(self$one_t(k = k)))
      }
      eic_all <- do.call(cbind, eic_all)
      return(eic_all)
    },
    clever_covariate = function(k) {
      if (self$A_intervene == 1) g <- self$g1W  else g <- 1 - self$g1W
      h_list <- list()
      for (t in 1:max(self$T_tilde)) {
        if (t > k) {
          # clever covariate is zero beyond
          h <- rep(0, length(g))
        } else {
          h <- -as.numeric(self$A == self$A_intervene) / g /
            self$density_censor$survival[, t] *
            self$density_failure$survival[, k] / self$density_failure$survival[, t]
        }
        h_list <- c(h_list, list(h))
      }
      # the first row is 1 ~ t_max for the first subject
      h_list <- do.call(cbind, h_list)
      # the first 1 ~ t_max element is for the first subject
      return(as.vector(t(h_list)))
    }
  )
)

# WILSON: what is G(t_ | xxx) ? I naively used survival function

#' compute multivariate normal quantile from a correlation matrix
#' @param corr correlation matrix
#' @param B number of monte-carlo samples drawn to estimate the quantile
#' @param alpha significant level (to compute 1-alpha quantile)
#' @return univariate numeric quantile of the quantile(max_j(abs(x))) where x is
#'  drawn from normal(0, corr)
#'
#' @export
compute_q <- function(corr, B = 1e3, alpha = 0.05) {
  dim <- nrow(corr)
  z <- apply(
    abs(MASS::mvrnorm(B, mu = rep(0, dim), Sigma = corr)), 1, max
  )
  return(as.numeric(stats::quantile(z, 1 - alpha)))
}
#' compute simutaneous confidence band around a survival curve esimator
#' @param eic_fit a matrix of efficient influence curve from `eic` class
#'  `all_t` method
#' @return a vector of standard error corresponding to each time point on the
#' survival curve
#'
#' @export
compute_simultaneous_ci <- function(eic_fit) {
  # compute the value to +- around the Psi_n
  n <- nrow(eic_fit)
  sigma_squared <- stats::cov(eic_fit)
  sigma <- stats::cor(eic_fit)
  # impute when the variance are zero
  sigma_squared[is.na(sigma_squared)] <- 1e-10
  sigma[is.na(sigma)] <- 1e-10

  variance_marginal <- diag(sigma_squared)
  q <- compute_q(corr = sigma, B = 1e3, alpha = 0.05)
  return(sqrt(variance_marginal) / sqrt(n) * q)
}

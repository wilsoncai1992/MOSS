require("R6")
require("SuperLearner")

#' onestep TMLE of treatment-rule specific survival curve
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods
#' @format \code{\link{R6Class}} object.
#' @examples
#' # MOSS_difference$new(dat, epsilon.step = 1e-5, max.iter = 1e3, tol = 1/nrow(dat), T.cutoff = NULL, verbose = FALSE)
#' @field dat data.frame
#' @field epsilon.step float
#' @field max.iter int
#' @field tol float
#' @field T.cutoff int
#' @field verbose bool
#' @section Methods:
#' @export
MOSS <- R6Class("MOSS",
  public = list(
    A = NULL,
    T_tilde = NULL,
    Delta = NULL,
    density_failure = NULL,
    density_censor = NULL,
    g1W = NULL,
    A_intervene = NULL,

    epsilon = NULL,
    max_num_interation = NULL,
    tmle_tolerance = NULL,
    k_grid = NULL,

    q_best = NULL,
    initialize = function(
      A,
      T_tilde,
      Delta,
      density_failure,
      density_censor,
      g1W,
      A_intervene,
      k_grid = NULL
    ) {
      self$A <- A
      self$T_tilde <- T_tilde
      self$Delta <- Delta
      self$density_failure <- density_failure
      self$density_censor <- density_censor
      self$g1W <- g1W
      self$A_intervene <- A_intervene

      self$k_grid <- k_grid
      return(self)
    },
    onestep_update_curve = function(density_failure, eic_fit, epsilon) {
      # don't handle pdf sum > 1
      pdf <- density_failure$pdf
      # simplify version 1
      # pdf2 <- pdf * exp(epsilon * eic_fit)

      # version 2: status quo
      # mean_eic <- colMeans(eic_fit)
      # v2 <- sqrt(sum(mean_eic ^ 2) / length(mean_eic))
      # # multiply - abs(mean_eic) to each row of the eic matrix
      # v1 <- t(- abs(mean_eic) * t(eic_fit))
      # pdf2 <- pdf * exp(epsilon * v1 / v2)

      # version 3: mark paper
      mean_eic <- colMeans(eic_fit)
      v2 <- sqrt(sum(mean_eic ^ 2))
      # multiply - abs(mean_eic) to each row of the eic matrix
      v1 <- t(- abs(mean_eic) * t(eic_fit))
      v1 <- apply(v1, 1, sum)
      pdf2 <- pdf * (1 + epsilon * v1 / v2) # more respecting LLFM

      density_failure2 <- survival_curve$new(t = density_failure$t, pdf = pdf2)
      density_failure2$pdf_to_survival()
      density_failure2$pdf_to_hazard()
      return(density_failure2)
    },
    onestep_curve = function(
      epsilon = 1e-5,
      max_num_interation = 1e2,
      tmle_tolerance = NULL,
      verbose = FALSE
    ) {
      self$epsilon <- epsilon
      self$max_num_interation <- max_num_interation
      if (is.null(tmle_tolerance)) {
        self$tmle_tolerance <- 1 / self$density_failure$n()
      } else {
        self$tmle_tolerance <- tmle_tolerance
      }

      k_grid <- 1:max(self$T_tilde)
      # create pdf
      self$density_failure <- self$density_failure$hazard_to_pdf()
      psi_n <- colMeans(self$density_failure$survival)
      eic_fit <- eic$new(
        A = self$A,
        T_tilde = self$T_tilde,
        Delta = self$Delta,
        density_failure = self$density_failure,
        density_censor = self$density_censor,
        g1W = self$g1W,
        psi = psi_n,
        A_intervene = self$A_intervene
      )$all_t(k_grid = k_grid)
      mean_eic <- colMeans(eic_fit)

      num_iteration <- 0

      mean_eic_inner_prod_prev <- abs(sqrt(sum(mean_eic ^ 2)))
      mean_eic_inner_prod_current <- mean_eic_inner_prod_prev
      mean_eic_inner_prod_best <- sqrt(sum(mean_eic ^ 2))
      self$q_best <- self$density_failure$clone(deep = TRUE)

      while (
        mean_eic_inner_prod_current >= self$tmle_tolerance * sqrt(max(k_grid))
      ) {
        if (verbose) {
          df_debug <- data.frame(num_iteration, mean_eic_inner_prod_current, mean(psi_n))
          colnames(df_debug) <- NULL
          print(df_debug)
        }
        # update
        self$density_failure <- self$onestep_update_curve(
          density_failure = self$density_failure,
          eic_fit = eic_fit,
          epsilon = self$epsilon
        )
        psi_n <- colMeans(self$density_failure$survival)
        eic_fit <- eic$new(
          A = self$A,
          T_tilde = self$T_tilde,
          Delta = self$Delta,
          density_failure = self$density_failure,
          density_censor = self$density_censor,
          g1W = self$g1W,
          psi = psi_n,
          A_intervene = self$A_intervene
        )$all_t(k_grid = k_grid)
        mean_eic <- colMeans(eic_fit)
        # new stopping
        mean_eic_inner_prod_prev <- mean_eic_inner_prod_current
        mean_eic_inner_prod_current <- abs(sqrt(sum(mean_eic ^ 2)))
        num_iteration <- num_iteration + 1

        # browser()
        if (mean_eic_inner_prod_prev < mean_eic_inner_prod_current) {
          self$epsilon <- - self$epsilon
        }
        if (mean_eic_inner_prod_current < mean_eic_inner_prod_best) {
          # the update caused PnEIC to beat the current best
          # update our best candidate
          self$q_best <- self$density_failure$clone(deep = TRUE)
          mean_eic_inner_prod_best <- mean_eic_inner_prod_current
        }

        if (num_iteration == self$max_num_interation) {
          break()
          warning("Max number of iteration reached, stop TMLE")
        }
      }
      # always output the best candidate for final result
      self$density_failure <- self$q_best
      psi_n <- colMeans(self$density_failure$survival)
      if (verbose) {
        message(paste(
          "Pn(EIC)=",
          formatC(mean_eic_inner_prod_best, format = "e", digits = 2),
          "Psi=",
          formatC(mean(psi_n), format = "e", digits = 2)
        ))
      }
      return(psi_n)
    },
    print = function() {
    }
  )
)


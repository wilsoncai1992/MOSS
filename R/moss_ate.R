#' onestep TMLE of average treatment effect on survival probabilities
#'
#' updating the hazard using constrained step size update
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data
#' @return Object of \code{\link{R6Class}} with methods
#' @format \code{\link{R6Class}} object.
#' @examples
#' # MOSS_hazard_ate$new(A = A, T_tilde = T.tilde, Delta = Delta, density_failure, density_censor, density_failure_0, density_censor_0, g1W, A_intervene = 1, k_grid = 1:max(T_tilde))
#' @field A vector of treatment
#' @field T_tilde vector of last follow up time
#' @field Delta vector of censoring indicator
#' @field density_failure survival_curve object of predicted counterfactual
#'  survival curve
#' @field density_censor survival_curve object of predicted counterfactual
#'  failure event survival curve
#' @field density_failure_0 survival_curve object of predicted counterfactual
#'  survival curve 0
#' @field density_censor_0 survival_curve object of predicted counterfactual
#'  failure event survival curve 0
#' @field g1W propensity score
#' @field k_grid vector of interested time points
#' @section Methods:
#' iterate_onestep update the initial estimator
#' @export
MOSS_hazard_ate <- R6Class("MOSS_hazard_ate",
  inherit = MOSS_hazard,
  public = list(
    density_failure_0 = NULL,
    density_censor_0 = NULL,

    q_best_0 = NULL,
    initialize = function(
      density_failure_0,
      density_censor_0,
      ...
    ) {
      super$initialize(...)
      self$density_failure_0 <- density_failure_0
      self$density_censor_0 <- density_censor_0
      return(self)
    },
    fit_epsilon = function(clipping = Inf) {
      dNt <- self$create_dNt()
      h_matrix_1 <- self$construct_long_data(
        A_intervene = 1,
        density_failure = self$density_failure,
        density_censor = self$density_censor
      )
      h_matrix_0 <- self$construct_long_data(
        A_intervene = 0,
        density_failure = self$density_failure_0,
        density_censor = self$density_censor_0
      )
      h_matrix <- h_matrix_1 - h_matrix_0

      offset_submodel_1 <- logit(as.vector(t(self$density_failure$hazard)))
      offset_submodel_0 <- logit(as.vector(t(self$density_failure_0$hazard)))
      offset_submodel <- offset_submodel_1
      offset_submodel[self$A == 0] <- offset_submodel_0[self$A == 0]
      submodel_fit <- glm.fit(
        x = h_matrix,
        y = dNt,
        family = binomial(),
        offset = offset_submodel,
        intercept = FALSE
      )
      epsilon_n <- submodel_fit$coefficients
      l2_norm <- sqrt(sum(epsilon_n ^ 2))
      if (l2_norm >= clipping) {
        # clipping the step size
        epsilon_n <- epsilon_n / l2_norm * clipping
      }

      hazard_new_1 <- expit(
        offset_submodel_1 +
        as.vector(h_matrix_1 %*% epsilon_n)
      )
      hazard_new_1 <- matrix(
        hazard_new_1,
        nrow = length(self$A),
        ncol = max(self$T_tilde),
        byrow = TRUE
      )
      hazard_new_0 <- expit(
        offset_submodel_0 +
        as.vector(h_matrix_0 %*% epsilon_n)
      )
      hazard_new_0 <- matrix(
        hazard_new_0,
        nrow = length(self$A),
        ncol = max(self$T_tilde),
        byrow = TRUE
      )
      # the new hazard for failure
      return(
        list(
          hazard_new_1 = survival_curve$new(
            t = 1:max(self$T_tilde), hazard = hazard_new_1
          )$hazard_to_survival(),
          hazard_new_0 = survival_curve$new(
            t = 1:max(self$T_tilde), hazard = hazard_new_0
          )$hazard_to_survival()
        )
      )
    },
    compute_mean_eic = function(psi_n, k_grid) {
      eic_fit_1 <- eic$new(
        A = self$A,
        T_tilde = self$T_tilde,
        Delta = self$Delta,
        density_failure = self$density_failure,
        density_censor = self$density_censor,
        g1W = self$g1W,
        psi = psi_n,
        A_intervene = 1
      )$all_t(k_grid = k_grid)
      eic_fit_0 <- eic$new(
        A = self$A,
        T_tilde = self$T_tilde,
        Delta = self$Delta,
        density_failure = self$density_failure_0,
        density_censor = self$density_censor_0,
        g1W = self$g1W,
        psi = psi_n,
        A_intervene = 0
      )$all_t(k_grid = k_grid)
      eic_fit <- eic_fit_1 - eic_fit_0
      mean_eic <- colMeans(eic_fit)
      return(mean_eic)
    },
    iterate_onestep = function(
      epsilon = 1e-2,
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

      psi_n <- colMeans(
        self$density_failure$survival - self$density_failure_0$survival
      )
      mean_eic <- self$compute_mean_eic(psi_n = psi_n, k_grid = k_grid)

      num_iteration <- 0

      mean_eic_inner_prod_prev <- abs(sqrt(sum(mean_eic ^ 2)))
      mean_eic_inner_prod_current <- mean_eic_inner_prod_prev
      mean_eic_inner_prod_best <- sqrt(sum(mean_eic ^ 2))
      self$q_best <- self$density_failure$clone(deep = TRUE)
      self$q_best_0 <- self$density_failure_0$clone(deep = TRUE)

      to_iterate <- TRUE
      if (is.infinite(mean_eic_inner_prod_current) | is.na(mean_eic_inner_prod_current)) {
        to_iterate <- FALSE
      }
      while (
        mean_eic_inner_prod_current >= self$tmle_tolerance * sqrt(max(k_grid)) &
        to_iterate
      ) {
        if (verbose) {
          df_debug <- data.frame(num_iteration, mean_eic_inner_prod_current, mean(psi_n))
          colnames(df_debug) <- NULL
          print(df_debug)
        }
        # update
        new_hazard <- self$fit_epsilon(clipping = self$epsilon)
        self$density_failure <- new_hazard$hazard_new_1
        self$density_failure_0 <- new_hazard$hazard_new_0

        psi_n <- colMeans(
          self$density_failure$survival - self$density_failure_0$survival
        )
        mean_eic <- self$compute_mean_eic(psi_n = psi_n, k_grid = k_grid)

        # new stopping
        mean_eic_inner_prod_prev <- mean_eic_inner_prod_current
        mean_eic_inner_prod_current <- abs(sqrt(sum(mean_eic ^ 2)))
        num_iteration <- num_iteration + 1
        if (is.infinite(mean_eic_inner_prod_current) | is.na(mean_eic_inner_prod_current)) {
          warning("stopping criteria diverged. Reporting best result so far.")
          break()
        }
        if (mean_eic_inner_prod_current < mean_eic_inner_prod_best) {
          # the update caused PnEIC to beat the current best
          # update our best candidate
          self$q_best <- self$density_failure$clone(deep = TRUE)
          self$q_best_0 <- self$density_failure_0$clone(deep = TRUE)
          mean_eic_inner_prod_best <- mean_eic_inner_prod_current
        }
        if (num_iteration == self$max_num_interation) {
          break()
          warning("Max number of iteration reached, stop TMLE")
        }
      }
      # always output the best candidate for final result
      self$density_failure <- self$q_best
      self$density_failure_0 <- self$q_best_0
      psi_n <- colMeans(
        self$density_failure$survival - self$density_failure_0$survival
      )
      if (verbose) {
        message(paste(
          "Pn(EIC)=",
          formatC(mean_eic_inner_prod_best, format = "e", digits = 2),
          "Psi=",
          formatC(mean(psi_n), format = "e", digits = 2)
        ))
      }
      return(psi_n)
    }
  )
)

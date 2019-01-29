require("R6")
require("SuperLearner")

#' onestep TMLE
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
MOSS_hazard <- R6Class("MOSS_hazard",
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
    create_dNt = function() {
      dNt <- matrix(0, nrow = length(self$A), ncol = max(self$T_tilde))
      for (i in 1:length(self$A)) {
        if (self$Delta[i] == 1) {
          dNt[i, self$T_tilde[i]] <- 1
        }
      }
      return(as.vector(t(dNt)))
    },
    construct_long_data = function() {
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
      )
      k_grid <- 1:max(self$T_tilde)
      h_matrix <- list()
      for (k in k_grid) {
        h <- eic_fit$clever_covariate(k = k)
        h_matrix <- c(h_matrix, list(h))
      }
      h_matrix <- do.call(cbind, h_matrix)
      return(h_matrix)
    },
    fit_epsilon = function(clipping = Inf) {
      dNt <- self$create_dNt()
      h_matrix <- self$construct_long_data()
      logit <- function(x) log(x) - log(1 - x)
      expit <- function(x) exp(x) / (1 + exp(x))
      submodel_fit <- glm.fit(
        x = h_matrix,
        y = dNt,
        family = binomial(),
        offset = logit(as.vector(t(self$density_failure$hazard))),
        intercept = FALSE
      )
      epsilon_n <- submodel_fit$coefficients
      l2_norm <- sqrt(sum(epsilon_n ^ 2))
      if (l2_norm >= clipping) {
        # clipping the step size
        epsilon_n <- epsilon_n / l2_norm * clipping
      }
      hazard_new <- expit(
        logit(as.vector(t(self$density_failure$hazard))) +
        as.vector(h_matrix %*% epsilon_n)
      )
      hazard_new <- matrix(
        hazard_new,
        nrow = length(self$A),
        ncol = max(self$T_tilde),
        byrow = TRUE
      )
      # the new hazard for failure
      return(
        survival_curve$new(
          t = 1:max(self$T_tilde), hazard = hazard_new
        )$hazard_to_survival()
      )
    },
    iterate_onestep = function(
      epsilon = 1e0,
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
        self$density_failure <- self$fit_epsilon(clipping = self$epsilon)

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
    }
  )
)


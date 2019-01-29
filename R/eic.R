require("R6")

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

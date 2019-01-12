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
      g <- ifelse(self$A_intervene == 1, self$g1W, 1 - self$g1W)
      part1_sum <- rep(0, length(g))
      for (t in 1:k) {
        h <- - as.numeric(self$A == self$A_intervene) / g /
          self$density_censor$survival[, t] *
          self$density_failure$survival[, k] / self$density_failure$survival[, t]
        part1 <-  h * (
          as.numeric(self$T_tilde == t & self$Delta == 1) -
          as.numeric(self$T_tilde >= t) * self$density_failure$hazard[, t]
        )
        part1_sum <- part1_sum + part1
      }
      part2 <- self$density_failure$survival[, k] - self$psi[, k]
      return(part1_sum + part2)
    },
    all_t = function() {

    }
  )
)

# WILSON: what is G(t_ | xxx) ? I naively used survival function
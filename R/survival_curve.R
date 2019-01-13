require("R6")

#' @export
survival_curve <- R6Class("survival_curve",
  public = list(
    t = NULL,
    hazard = NULL,
    # P( T >= t)
    survival = NULL,
    initialize = function(t, hazard = NULL, survival = NULL) {
      # only supports integer grid
      from_hazard <- !is.null(hazard)
      from_survival <- !is.null(survival)
      if (from_hazard & from_survival) stop("cannot construct from both")
      if (!all.equal(t, seq(range(t)[1], range(t)[2]))) {
        stop("t is not integer without gap")
      }
      self$t <- t
      if (from_hazard) {
        message("construct from hazard")
        self$hazard <- as.matrix(hazard)
      }
      if (from_survival) {
        message("construct from survival")
        self$survival <- as.matrix(survival)
      }
    },
    n = function() {
      n1 <- nrow(self$hazard)
      n2 <- nrow(self$survival)
      return(ifelse(is.null(n1), n2, n1))
    },
    hazard_to_survival = function() {
      self$survival <- matrix(NA, nrow = self$n(), ncol = max(self$t))
      # for (i in 1:self$n()) {
      #   self$survival[i, ] <- cumprod(1 - self$hazard[i, ])
      # }

      # cumulative hazard approach
      # integral from t = 0 and not include right bound. e.g. [0, 1)
      hazard_integral <- cbind(0, self$hazard)
      hazard_integral <- hazard_integral[, -ncol(hazard_integral)]
      for (i in 1:self$n()) {
        self$survival[i, ] <- exp(- cumsum(hazard_integral[i, ]))
      }
    },
    pdf_to_survival = function() {

    },
    pdf_survival_to_hazard = function() {

    },
    display = function() {

    }
  )
)

# S ~ A = 1, W, t plot



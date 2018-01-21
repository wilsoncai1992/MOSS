require("R6")
require("SuperLearner")

#' @export
MOSS_difference <- R6Class("MOSS_difference",
  public = list(
    # simultaneous CI
    simul_CI = NULL,
    initialize = function(dat,
                          dW,
                          epsilon.step = 1e-5,
                          max.iter = 1e3,
                          tol = 1/nrow(dat),
                          T.cutoff = NULL,
                          verbose = FALSE) {
      self$dat <- dat
      self$dW <- dW
      self$epsilon.step <- epsilon.step
      self$max.iter <- max.iter
      self$tol <- tol
      self$T.cutoff <- T.cutoff
      self$verbose <- verbose

      self$check_and_preprocess_data(T.cutoff = self$T.cutoff)
      self$update_tensor <- matrix(0, nrow = self$n_sample, ncol = length(self$T.uniq))
      # self$inside_exp <- rep(0, length(self$T.uniq))
      self$inside_exp <- matrix(0, ncol = length(self$T.uniq), nrow = self$n_sample)

      self$MOSS_A1 <- onestepfit = MOSS$new(dat, dW = 1,
                verbose = self$verbose, epsilon.step = epsilon.step, max.iter = max.iter)
      self$MOSS_A0 <- onestepfit = MOSS$new(dat, dW = 0,
                verbose = self$verbose, epsilon.step = epsilon.step, max.iter = max.iter)

    },
    initial_fit = function(g.SL.Lib = c("SL.mean","SL.glm", 'SL.gam'),
                           Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam"),
                           ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam")){
      self$MOSS_A1$initial_fit(g.SL.Lib = g.SL.Lib,
                               Delta.SL.Lib = Delta.SL.Lib,
                               ht.SL.Lib = ht.SL.Lib)
      self$MOSS_A1$transform_failure_hazard_to_survival()
      self$MOSS_A1$transform_failure_hazard_to_pdf()
      self$MOSS_A1$compute_EIC()
      self$MOSS_A0$initial_fit(g.SL.Lib = g.SL.Lib,
                               Delta.SL.Lib = Delta.SL.Lib,
                               ht.SL.Lib = ht.SL.Lib)
      self$MOSS_A0$transform_failure_hazard_to_survival()
      self$MOSS_A0$transform_failure_hazard_to_pdf()
      self$MOSS_A0$compute_EIC()
    }
  )
)

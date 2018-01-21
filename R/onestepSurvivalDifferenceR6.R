require("R6")
require("SuperLearner")

#' @export
MOSS_difference <- R6Class("MOSS_difference",
  public = list(
    dat = NULL,
    epsilon.step = NULL,
    max.iter = NULL,
    tol = NULL,
    T.cutoff = NULL,
    verbose = NULL,
# single MOSS object
    MOSS_A1 = NULL,
    MOSS_A0 = NULL,
# targeting
    D1.t = NULL,
    Pn.D1.t = NULL,
    stopping_history = numeric(),
    Psi.hat = NULL,
    # simultaneous CI
    simul_CI = NULL,
    initialize = function(dat,
                          epsilon.step = 1e-5,
                          max.iter = 1e3,
                          tol = 1/nrow(dat),
                          T.cutoff = NULL,
                          verbose = FALSE) {

      self$MOSS_A1 <- MOSS$new(dat, dW = 1,
                              epsilon.step = epsilon.step,
                              max.iter = max.iter,
                              tol = tol,
                              T.cutoff = T.cutoff,
                              verbose = verbose)
      self$MOSS_A0 <- MOSS$new(dat, dW = 0,
                              epsilon.step = epsilon.step,
                              max.iter = max.iter,
                              tol = tol,
                              T.cutoff = T.cutoff,
                              verbose = verbose)

      self$dat <- self$MOSS_A1$dat
      self$epsilon.step <- self$MOSS_A1$epsilon.step
      self$max.iter <- self$MOSS_A1$max.iter
      self$tol <- self$MOSS_A1$tol
      self$T.cutoff <- self$MOSS_A1$T.cutoff
      self$verbose <- self$MOSS_A1$verbose
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
    },
    onestep_diff_curve = function(){
      browser()
      self$D1.t <- self$MOSS_A1$D1.t - self$MOSS_A0$D1.t
      self$Pn.D1.t <- colMeans(self$D1.t)
      self$MOSS_A1$D1.t <- self$D1.t
      self$MOSS_A0$D1.t <- self$D1.t
      self$MOSS_A1$Pn.D1.t <- self$Pn.D1.t
      self$MOSS_A0$Pn.D1.t <- self$Pn.D1.t

      iter_count <- 0
      stopping_prev <- Inf

      stopping <- self$MOSS_A1$compute_stopping()
      while ((stopping >= self$tol) & (iter_count <= self$max.iter)) {
      # while ((stopping >= self$tol) & (iter_count <= self$max.iter) & ((stopping_prev - stopping) >= max(-self$tol, -1e-5))) {
        print(stopping)
        if (stopping_prev < stopping) onestepfit$epsilon.step <- -onestepfit$epsilon.step
        # self$onestep_curve_update()
        self$MOSS_A1$onestep_curve_update_mat()
        self$MOSS_A0$onestep_curve_update_mat()
        self$MOSS_A1$compute_EIC()
        self$MOSS_A0$compute_EIC()

        self$D1.t <- self$MOSS_A1$D1.t - self$MOSS_A0$D1.t
        self$Pn.D1.t <- colMeans(self$D1.t)
        self$MOSS_A1$D1.t <- self$D1.t
        self$MOSS_A0$D1.t <- self$D1.t
        self$MOSS_A1$Pn.D1.t <- self$Pn.D1.t
        self$MOSS_A0$Pn.D1.t <- self$Pn.D1.t

        iter_count <- iter_count + 1
        self$stopping_history[iter_count] <- stopping
        stopping_prev <- self$stopping_history[iter_count]
        stopping <- self$MOSS_A1$compute_stopping()

        # if (iter_count %% 10 == 0) onestepfit$plot_onestep_curve(add = TRUE)
        # if (iter_count %% 10 == 0) plot(onestepfit$Pn.D1.t); abline(h = 0)
      }

      if (iter_count == self$max.iter) {
        warning('Max Iter count reached, stop iteration.')
      }

      self$MOSS_A1$compute_Psi()
      self$MOSS_A0$compute_Psi()
      self$Psi.hat <- self$MOSS_A1$Psi.hat - self$MOSS_A0$Psi.hat
    },
    print = function(){
      out <- data.frame(self$T.uniq, self$Psi.hat, self$sd_EIC, self$upper_CI, self$lower_CI)
      colnames(out) <- c('Time', 'survival curve', 'std_err', 'upper_CI', 'lower_CI')
      return(out)
    },
  )
)

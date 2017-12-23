require("R6")
require("SuperLearner")

#' @export
MOSS <- R6Class("MOSS",
  public = list(
    dat = NULL,
    dW = NULL,
    epsilon.step = NULL,
    max.iter = NULL,
    tol = NULL,
    T.cutoff = NULL,
    verbose = NULL,

    n_sample = NULL,
    W_names = NULL,
    W = NULL,
    A = NULL,
    T.tilde = NULL,
    Delta = NULL,
    T.uniq = NULL,
    K = NULL,
    T.max = NULL,

    # g
    g.fitted = NULL,
    # N(t)
    h.hat.t = NULL,
    h.hat.t_full = NULL,
    Qn.A1.t_full = NULL,
    Qn.A1.t = NULL,
    qn.A1.t_full = NULL,
    qn.A1.t = NULL,
    # A_c(t)
    Gn.A1.t_full = NULL,
    Gn.A1.t = NULL,
    # EIC
    D1.t = NULL,
    D1.A1.t = NULL,
    Pn.D1.t = NULL,
    # targeting
    stopping_criteria = NULL,
    stopping_history = numeric(),
    update_tensor = NULL,
    inside_exp = 0,
    Psi.hat = NULL,
    sd_EIC = NULL,
    upper_CI = NULL,
    lower_CI = NULL,
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
    },
    check_and_preprocess_data = function(nbin = 4, T.cutoff = NULL){
      message('check data validity')
      self$W_names <- grep('W', colnames(self$dat), value = TRUE)
      # colnames should exist
      if (!('T.TILDE' %in% toupper(colnames(self$dat)))) stop('T.tilde should exist')
      if (!('A' %in% colnames(self$dat))) stop('A should exist')
      if (!('Delta' %in% colnames(self$dat))) {
        warning('Delta not found. Set Delta = 1.')
        self$dat$Delta <- rep(1, nrow(self$dat))
      }
      # keep necessary columns
      self$dat <- self$dat[,c('T.tilde', 'A', 'Delta', self$W_names)]
      # dW length should be same
      if (is.null(self$dW)) stop('Input dW!')
      if (length(self$dW) == 1) self$dW <- rep(self$dW, nrow(self$dat))
      if (length(self$dW) != nrow(self$dat)) stop('Input dW should have same length as dat!')

      # keep only T.tilde > 0
      to_keep <- self$dat$T.tilde != 0
      self$dW <- self$dW[to_keep]
      self$dat <- self$dat[to_keep,]

      self$n_sample <- nrow(self$dat)
      # number of samples should be same
      if (length(self$dW) != self$n_sample) stop('The length of input dW is not same as the sample size!')

      # create objects
      self$W <- self$dat[,self$W_names]
      self$W <- as.data.frame(self$W)
      self$A <- self$dat$A
      self$Delta <- self$dat$Delta
      self$T.tilde <- self$dat$T.tilde
      self$T.uniq <- sort(unique(self$dat$T.tilde))
      self$K <- length(self$T.uniq)
      self$T.max <- max(self$T.uniq)
    },
    initial_fit = function(g.SL.Lib = c("SL.mean","SL.glm", "SL.step", "SL.glm.interaction"),
                           Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                           ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth")){
      # browser()
      message('initial fit')
      fit_out <- initial_SL_fit(ftime = self$T.tilde,
                                ftype = self$Delta,
                                trt = self$A,
                                adjustVars = data.frame(self$W),
                                t_0 = self$T.max,
                                SL.trt = g.SL.Lib,
                                SL.ctime = Delta.SL.Lib,
                                SL.ftime = ht.SL.Lib)
      haz1 <- fit_out[[1]]
      haz0 <- fit_out[[2]]
      S_Ac_1 <- fit_out[[3]]
      S_Ac_0 <- fit_out[[4]]
      g_1 <- fit_out[[5]]
      g_0 <- fit_out[[6]]

      if (all(self$dW == 1)) haz <- haz1; S_Ac <- S_Ac_1; self$g.fitted <- g_1
      if (all(self$dW == 0)) haz <- haz0; S_Ac <- S_Ac_0; self$g.fitted <- g_0

      self$h.hat.t_full <- as.matrix(haz)
      self$h.hat.t <- self$h.hat.t_full[,self$T.uniq]
      self$Gn.A1.t_full <- as.matrix(S_Ac)
      self$Gn.A1.t <- self$Gn.A1.t_full[,self$T.uniq]
    },
    transform_failure_hazard_to_survival = function(){
      Qn.A1.t <- matrix(0, nrow = self$n_sample, ncol = length(self$T.uniq))
      Qn.A1.t_full <- matrix(NA, nrow = self$n_sample, ncol = ncol(self$h.hat.t_full))
      # cum-product approach (2016-10-05)
      for (it in 1:self$n_sample) {
        Qn.A1.t_full[it,] <- cumprod(1 - self$h.hat.t_full[it,])
      }
      self$Qn.A1.t_full <- Qn.A1.t_full
      self$Qn.A1.t <- Qn.A1.t_full[, self$T.uniq]
    },
    transform_failure_hazard_to_pdf = function(){
      qn.A1.t_full <- matrix(0, nrow = self$n_sample, ncol = ncol(self$Qn.A1.t_full))
      for (it in 1:self$n_sample) {
        qn.A1.t_full[it,] <- self$h.hat.t_full[it,] * self$Qn.A1.t_full[it,]
      }
      self$qn.A1.t_full <- qn.A1.t_full
      self$qn.A1.t <- qn.A1.t_full[, self$T.uniq]

      # self$qn.A1.t_full <- survivalDensity$new(pdf = discreteDensity$new(p = qn.A1.t_full, t_grid = self$T.uniq))
      # self$qn.A1.t <- survivalDensity$new(pdf = discreteDensity$new(p = qn.A1.t_full[, self$T.uniq], t_grid = self$T.uniq))

    },
    compute_EIC = function(){
      # browser()
      I.A.dW <- self$A == self$dW

      # D_1* in paper
      D1.t <- matrix(0, nrow = self$n_sample, ncol = self$K)

      for (it in 1:self$n_sample) {
        t_Delta1.vec <- create_Yt_vector_with_censor(Time = self$T.tilde[it],
                                                     Delta = self$Delta[it],
                                                     t.vec = 1:self$T.max)
        t.vec <- create_Yt_vector(Time = self$T.tilde[it],
                                  t.vec = 1:self$T.max)
        alpha2 <- t_Delta1.vec - t.vec * self$h.hat.t_full[it,]
        h1 <- -I.A.dW[it]/self$g.fitted[it]/self$Gn.A1.t_full[it,]/self$Qn.A1.t_full[it,]

        not_complete <- h1 * alpha2
        # D1 matrix
        D1.t[it, ] <- cumsum(not_complete)[self$T.uniq] * self$Qn.A1.t_full[it,self$T.uniq] # complete influence curve
      }

      # turn unstable results to 0
      D1.t[is.na(D1.t)] <- 0

      self$D1.t <- D1.t
      self$Pn.D1.t <- colMeans(self$D1.t)
    },
    compute_stopping = function(){
      return(sqrt(l2_inner_prod_step(self$Pn.D1.t, self$Pn.D1.t, self$T.uniq))/length(self$T.uniq))
    },
    compute_hazard_from_pdf_and_survival = function(){
      hazard_new <- matrix(0, nrow = self$n_sample, ncol = self$T.max)
      for (it in 1:self$n_sample) {
        hazard_new[it, ] <- self$qn.A1.t_full[it, ] / self$Qn.A1.t_full[it,]
      }
      # dirty fix: upper bound hazard
      hazard_new[hazard_new >= 1] <- .8

      self$h.hat.t_full <- hazard_new
      self$h.hat.t <- hazard_new[, self$T.uniq]
    },
    compute_survival_from_pdf = function(){
      self$Qn.A1.t <- compute_step_cdf(pdf.mat = self$qn.A1.t, t.vec = self$T.uniq, start = Inf)
      self$Qn.A1.t_full <- compute_step_cdf(pdf.mat = self$qn.A1.t_full, t.vec = 1:self$T.max, start = Inf)
    },
    onestep_curve_update_pooled = function(){
      update <- compute_onestep_update_matrix(D1.t.func.prev = self$D1.t,
                                              Pn.D1.func.prev = self$Pn.D1.t,
                                              dat = self$dat,
                                              T.uniq = self$T.uniq,
                                              W_names = self$W_names,
                                              dW = self$dW)
      self$inside_exp <- sum(update)
      # self$inside_exp[is.na(self$inside_exp)] <- 0
      self$qn.A1.t <- self$qn.A1.t * exp(self$epsilon.step * self$inside_exp)
      self$qn.A1.t_full <- self$qn.A1.t_full * exp(self$epsilon.step * self$inside_exp)

      # For density sum > 1: normalize the updated qn
      # norm.factor <- compute_step_cdf(pdf.mat = self$qn.A1.t, t.vec = self$T.uniq, start = Inf)[,1]
      # self$qn.A1.t[norm.factor > 1,] <- self$qn.A1.t[norm.factor > 1,] / norm.factor[norm.factor > 1]
      # self$qn.A1.t_full[norm.factor > 1,] <- self$qn.A1.t_full[norm.factor > 1,] / norm.factor[norm.factor > 1]

      # # if some qn becomes all zero, prevent NA exisitence
      # self$qn.A1.t[is.na(self$qn.A1.t)] <- 0
      # self$qn.A1.t_full[is.na(self$qn.A1.t_full)] <- 0

      self$qn.A1.t_full <- self$qn.A1.t_full / rowSums(self$qn.A1.t_full)
      self$qn.A1.t <- self$qn.A1.t_full[,self$T.uniq]

      # compute new Survival
      self$compute_survival_from_pdf()

      # compute new hazard
      # self$compute_hazard_from_pdf_and_survival()
    },
    onestep_curve_update_mat = function(){
      update <- compute_onestep_update_matrix(D1.t.func.prev = self$D1.t,
                                              Pn.D1.func.prev = self$Pn.D1.t,
                                              dat = self$dat,
                                              T.uniq = self$T.uniq,
                                              W_names = self$W_names,
                                              dW = self$dW)
      self$inside_exp <- (update)

      self$qn.A1.t <- self$qn.A1.t * exp(self$epsilon.step * self$inside_exp)

      inside_exp_longer <- matrix(NA, ncol = self$T.max, nrow = self$n_sample)
      inside_exp_longer[,self$T.uniq] <- self$inside_exp
      inside_exp_longer <- t(zoo::na.locf(t(inside_exp_longer)))
      self$qn.A1.t_full <- self$qn.A1.t_full * exp(self$epsilon.step * inside_exp_longer)

      self$qn.A1.t_full <- self$qn.A1.t_full / rowSums(self$qn.A1.t_full)
      self$qn.A1.t <- self$qn.A1.t_full[,self$T.uniq]

      # compute new Survival
      self$compute_survival_from_pdf()

      # compute new hazard
      self$compute_hazard_from_pdf_and_survival()
    },
    onestep_curve_update = function(){
      update <- compute_onestep_update_matrix(D1.t.func.prev = self$D1.t,
                                              Pn.D1.func.prev = self$Pn.D1.t,
                                              dat = self$dat,
                                              T.uniq = self$T.uniq,
                                              W_names = self$W_names,
                                              dW = self$dW)
      self$inside_exp <- colSums(update)
      # self$inside_exp[is.na(self$inside_exp)] <- 0

      self$qn.A1.t <- multiple_vector_to_matrix(self$qn.A1.t, exp(self$epsilon.step * self$inside_exp))

      inside_exp_longer <- rep(NA, self$T.max)
      inside_exp_longer[self$T.uniq] <- self$inside_exp
      inside_exp_longer <- zoo::na.locf(inside_exp_longer)
      self$qn.A1.t_full <- multiple_vector_to_matrix(self$qn.A1.t_full, exp(self$epsilon.step * inside_exp_longer))

      self$qn.A1.t_full <- self$qn.A1.t_full / rowSums(self$qn.A1.t_full)
      self$qn.A1.t <- self$qn.A1.t_full[,self$T.uniq]

      # compute new Survival
      self$compute_survival_from_pdf()

      # compute new hazard
      self$compute_hazard_from_pdf_and_survival()
    },
    onestep_curve_update_no_normalize = function(){
      update <- compute_onestep_update_matrix(D1.t.func.prev = self$D1.t,
                                              Pn.D1.func.prev = self$Pn.D1.t,
                                              dat = self$dat,
                                              T.uniq = self$T.uniq,
                                              W_names = self$W_names,
                                              dW = self$dW)
      self$inside_exp <- colSums(update)
      # self$inside_exp[is.na(self$inside_exp)] <- 0

      self$qn.A1.t <- multiple_vector_to_matrix(self$qn.A1.t, exp(self$epsilon.step * self$inside_exp))

      inside_exp_longer <- rep(NA, self$T.max)
      inside_exp_longer[self$T.uniq] <- self$inside_exp
      inside_exp_longer <- zoo::na.locf(inside_exp_longer)
      self$qn.A1.t_full <- multiple_vector_to_matrix(self$qn.A1.t_full, exp(self$epsilon.step * inside_exp_longer))

      # self$qn.A1.t_full <- self$qn.A1.t_full / rowSums(self$qn.A1.t_full)
      self$qn.A1.t <- self$qn.A1.t_full[,self$T.uniq]

      # compute new Survival
      self$compute_survival_from_pdf()

      # compute new hazard
      self$compute_hazard_from_pdf_and_survival()
    },
    compute_Psi = function(){
      self$Psi.hat <- colMeans(self$Qn.A1.t_full)

      self$sd_EIC <- rep(NA, self$T.max)
      self$sd_EIC[self$T.uniq] <- sqrt(apply(self$D1.t, 2, var)/self$n_sample)
      self$sd_EIC <- zoo::na.locf(self$sd_EIC)
      self$upper_CI <- self$Psi.hat + 1.96 * self$sd_EIC
      self$lower_CI <- self$Psi.hat - 1.96 * self$sd_EIC

      EIC_sup_norm <- abs(self$Pn.D1.t)
    },
    onestep_curve = function(g.SL.Lib = c("SL.mean","SL.glm", 'SL.gam'),
                             Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam"),
                             ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam")){
      self$initial_fit(g.SL.Lib = g.SL.Lib,
                       Delta.SL.Lib = Delta.SL.Lib,
                       ht.SL.Lib = ht.SL.Lib)
      self$transform_failure_hazard_to_survival()
      self$transform_failure_hazard_to_pdf()
      self$compute_EIC()

      iter_count <- 0
      stopping_prev <- Inf
      all_loglikeli <- numeric()

      stopping <- self$compute_stopping()
      while ((stopping >= self$tol) & (iter_count <= self$max.iter)) {
      # while ((stopping >= self$tol) & (iter_count <= self$max.iter) & ((stopping_prev - stopping) >= max(-self$tol, -1e-5))) {
        print(stopping)
        if (stopping_prev < stopping) onestepfit$epsilon.step <- -onestepfit$epsilon.step
        # self$onestep_curve_update()
        self$onestep_curve_update_mat()
        self$compute_EIC()
        iter_count <- iter_count + 1
        self$stopping_history[iter_count] <- stopping
        stopping_prev <- self$stopping_history[iter_count]
        stopping <- self$compute_stopping()

        # if (iter_count %% 10 == 0) onestepfit$print_onestep_curve(add = TRUE)
        # if (iter_count %% 10 == 0) plot(onestepfit$Pn.D1.t); abline(h = 0)
      }

      if (iter_count == self$max.iter) {
        warning('Max Iter count reached, stop iteration.')
      }

      self$compute_Psi()
    },
    print_onestep_curve = function(...){
      step_curve <- stepfun(x = 1:self$T.max, y = c(self$Psi.hat, self$Psi.hat[length(self$Psi.hat)]))
      # can `add`, `col`
      curve(step_curve, from = 0, to = self$T.max, ...)
    },
    print = function(){
      data.frame(self$Psi.hat, self$sd_EIC, self$upper_CI, self$lower_CI)
    },
    plot_CI_pointwise = function(...){
      polygon(c(1:self$T.max, rev(1:self$T.max)), c(c(self$upper_CI), rev(c(self$lower_CI))),
                          col = rgb(0.7,0.7,0.7,0.4),
                          border = NA,
                          ...)
      self$print_onestep_curve(...)
    }
  )
)


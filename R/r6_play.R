require("R6")

#' @export
MOSS <- R6Class("MOSS",
  public = list(
    dat = NULL,
    dW = NULL,
    g.SL.Lib = NULL,
    Delta.SL.Lib = NULL,
    ht.SL.Lib = NULL,
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
    T.uniq = NULL,
    K = NULL,
    T.max = NULL,

    # g
    gHatSL = NULL,
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
    initialize = function(dat,
                          dW,
                          g.SL.Lib = c("SL.glm", "SL.step", "SL.glm.interaction"),
                          Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                          ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                          epsilon.step = 1e-5,
                          max.iter = 1e3,
                          tol = 1/nrow(dat),
                          T.cutoff = NULL,
                          verbose = FALSE) {
      self$dat <- dat
      self$dW <- dW
      self$g.SL.Lib <- g.SL.Lib
      self$Delta.SL.Lib <- Delta.SL.Lib
      self$ht.SL.Lib <- ht.SL.Lib
      self$epsilon.step <- epsilon.step
      self$max.iter <- max.iter
      self$tol <- tol
      self$T.cutoff <- T.cutoff
      self$verbose <- verbose

      self$check_and_preprocess_data()
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
        self$T.tilde <- self$dat$T.tilde
        self$T.uniq <- sort(unique(self$dat$T.tilde))
        self$K <- length(self$T.uniq)
        self$T.max <- max(self$T.uniq)

        if(all(self$dW == 0)) {
          self$A <- 1 - self$A
          self$dW <- 1 - self$dW
        }else if(all(self$dW == 1)){
        }else{
          stop('not implemented!')
        }
    },
    fit_g_initial = function(){
      message('fit g')
      self$gHatSL <- SuperLearner(Y = self$A,
                                 X = self$W,
                                 SL.library = self$g.SL.Lib,
                                 family = "binomial")
      # g.hat for each observation
      self$g.fitted <- self$gHatSL$SL.predict
    },
    fit_failure_hazard = function(){
      message('fit failure hazard')
      h.hat.t <- estimate_hazard_SL(dat = self$dat,
                                   T.uniq = self$T.uniq,
                                   ht.SL.Lib = self$ht.SL.Lib)
      # h.hat at all time t=[0,t.max]
      self$h.hat.t_full <- as.matrix(h.hat.t$out_haz_full)
      # h.hat at observed unique time t = T.grid
      self$h.hat.t <- as.matrix(h.hat.t$out_haz)
    },
    fit_censoring_cdf = function(cutoff = 0.05){
      message('fit censoring cdf')
      G.hat.t <- estimate_censoring_SL(dat = self$dat,
                                       T.uniq = self$T.uniq,
                                       Delta.SL.Lib = self$Delta.SL.Lib)
      if(any(G.hat.t$out_censor_full <= cutoff)){
        warning('G.hat has extreme small values! lower truncate to 0.05')
        G.hat.t$out_censor_full[G.hat.t$out_censor_full < cutoff] <- cutoff
        G.hat.t$out_censor[G.hat.t$out_censor < cutoff] <- cutoff
      }

      self$Gn.A1.t_full <- as.matrix(G.hat.t$out_censor_full)
      self$Gn.A1.t <- as.matrix(G.hat.t$out_censor)
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
    },
    compute_EIC = function(){
      I.A.dW <- self$A == self$dW

      D1.t <- matrix(0, nrow = self$n_sample, ncol = self$K)
      D1.A1.t <- matrix(0, nrow = self$n_sample, ncol = self$K)

      for (it in 1:self$n_sample) {
        t_Delta1.vec <- create_Yt_vector_with_censor(Time = self$T.tilde[it],
                                                     Delta = self$Delta[it],
                                                     t.vec = 1:self$T.max)
        t.vec <- create_Yt_vector(Time = self$T.tilde[it],
                                  t.vec = 1:self$T.max)
        alpha2 <- t_Delta1.vec - t.vec * self$h.hat.t_full[it,]

        alpha1 <- -I.A.dW[it]/self$g.fitted[it]/self$Gn.A1.t_full[it,]/self$Qn.A1.t_full[it,]
        alpha1_A1 <- -1/self$g.fitted[it]/self$Gn.A1.t_full[it,]/self$Qn.A1.t_full[it,]

        not_complete <- alpha1 * alpha2
        not_complete_A1 <- alpha1_A1 * alpha2
        # D1 matrix
        D1.t[it, ] <- cumsum(not_complete)[self$T.uniq] * self$Qn.A1.t[it,] # complete influence curve
        D1.A1.t[it, ] <- cumsum(not_complete_A1)[self$T.uniq] * self$Qn.A1.t[it,] # also update those A = 0.
      }

      # turn unstable results to 0
      D1.t[is.na(D1.t)] <- 0
      D1.A1.t[is.na(D1.A1.t)] <- 0
      self$D1.t <- D1.t
      self$D1.A1.t <- D1.A1.t

      self$Pn.D1.t <- colMeans(self$D1.t)
    },
    display = function() {
      cat("Make = ", self$make,
          " Price = ", self$price, "\n")
    }
  )
)


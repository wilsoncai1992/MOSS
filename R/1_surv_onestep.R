#' One-step TMLE estimator for survival curve
#'
#' one-step TMLE estimate of the treatment specific survival curve. Under right-censored data
#'
#' options to ADD:
#' SL.formula: the covariates to include in SL
#'
#' @param dat A data.frame with columns T.tilde, delta, A, W. T.tilde = min(T, C) is either the failure time of censor time, whichever happens first. 'delta'= I(T <= C) is the indicator of whether we observe failure time. A is binary treatment. W is baseline covariates. All columns with character "W" will be treated as baseline covariates.
#' @param dW A binary vector specifying dynamic treatment (as a function output of W)
#' @param g.SL.Lib A vector of string. SuperLearner library for fitting treatment regression
#' @param Delta.SL.Lib A vector of string. SuperLearner library for fitting censoring regression
#' @param ht.SL.Lib A vector of string. SuperLearner library for fitting conditional hazard regression
#' @param epsilon.step numeric. step size for one-step recursion
#' @param max.iter integer. maximal number of recursion for one-step
#' @param tol numeric. tolerance for optimization
#' @param T.cutoff int. Enforce randomized right-censoring to the observed data, so that don't estimate survival curve beyond a time point. Useful when time horizon is long.
#' @param verbose boolean. When TRUE, plot the initial fit curve, and output the objective function value during optimzation
#' @param ... additional options for plotting initial fit curve
#'
#' @return Psi.hat A numeric vector of estimated treatment-specific survival curve
#' @return T.uniq A vector of descrete time points where Psi.hat take values (have same length as Psi.hat)
#' @return params A list of estimation parameters set by user
#' @return variables A list of data summary statistics
#' @return initial_fit A list of initial fit values (hazard, g_1, Delta)
#'
#' @export
#'
#' @examples
#' library(simcausal)
#' D <- DAG.empty()
#'
#' D <- D +
#'     node("W", distr = "rbinom", size = 1, prob = .5) +
#'     node("A", distr = "rbinom", size = 1, prob = .15 + .5*W) +
#'     node("Trexp", distr = "rexp", rate = 1 + .5*W - .5*A) +
#'     node("Cweib", distr = "rweibull", shape = .7 - .2*W, scale = 1) +
#'     node("T", distr = "rconst", const = round(Trexp*100,0)) +
#'     node("C", distr = "rconst", const = round(Cweib*100, 0)) +
#'     node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
#'     node("delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
#' setD <- set.DAG(D)
#'
#' dat <- sim(setD, n=3e2)
#'
#' library(dplyr)
#' # only grab ID, W's, A, T.tilde, Delta
#' Wname <- grep('W', colnames(dat), value = TRUE)
#' dat <- dat[,c('ID', Wname, 'A', "T.tilde", "delta")]
#'
#' dW <- rep(1, nrow(dat))
#' onestepfit <- surv_onestep(dat = dat,
#'                             dW = dW,
#'                             verbose = FALSE,
#'                             epsilon.step = 1e-3,
#'                             max.iter = 1e3)
#' @import dplyr
#' @import survtmle
#' @import abind
#' @import SuperLearner
surv_onestep <- function(dat,
                         dW = rep(1, nrow(dat)),
                         g.SL.Lib = c("SL.glm", "SL.step", "SL.glm.interaction"),
                         Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                         ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                         epsilon.step = 1e-5,
                         max.iter = 1e3,
                         tol = 1/nrow(dat),
                         T.cutoff = NULL,
                         verbose = FALSE,
                         ...) {
  # ===================================================================================
  # preparation
  # ===================================================================================
  after_check <- check_and_preprocess_data(dat = dat, dW = dW, T.cutoff = T.cutoff)
  dat <- after_check$dat
  dW <- after_check$dW
  n.data <- after_check$n.data
  W_names <- after_check$W_names

  W <- dat[,W_names]
  W <- as.data.frame(W)

  # dW check
  if(all(dW == 0)) {
    dat$A <- 1 - dat$A # when dW is all zero
    dW <- 1 - dW
  }else if(all(dW == 1)){

  }else{
    stop('not implemented!')
  }

  T.uniq <- sort(unique(dat$T.tilde))
  T.max <- max(T.uniq)
  # ===================================================================================
  # estimate g(A|W)
  # ===================================================================================
  gHatSL <- SuperLearner(Y=dat$A, X=W, SL.library=g.SL.Lib, family="binomial")
  # g.hat for each observation
  g.fitted <- gHatSL$SL.predict
  # ===================================================================================
  # conditional hazard (by SL)
  # ===================================================================================
  message('estimating conditional hazard')

  h.hat.t <- estimate_hazard_SL(dat = dat, T.uniq = T.uniq, ht.SL.Lib = ht.SL.Lib)
  # h.hat at all time t=[0,t.max]
  h.hat.t_full <- as.matrix(h.hat.t$out_haz_full)
  # h.hat at observed unique time t = T.grid
  h.hat.t <- as.matrix(h.hat.t$out_haz)
  # ===================================================================================
  # estimate censoring G(A|W)
  # ===================================================================================
  message('estimating censoring')
  G.hat.t <- estimate_censoring_SL(dat = dat, T.uniq = T.uniq,
                                   Delta.SL.Lib = Delta.SL.Lib)
  # cutoff <- 0.1
  cutoff <- 0.05
  if(any(G.hat.t$out_censor_full <= cutoff)){
    warning('G.hat has extreme small values! lower truncate to 0.05')
    G.hat.t$out_censor_full[G.hat.t$out_censor_full < cutoff] <- cutoff
    G.hat.t$out_censor[G.hat.t$out_censor < cutoff] <- cutoff
  }

  Gn.A1.t_full <- as.matrix(G.hat.t$out_censor_full)
  Gn.A1.t <- as.matrix(G.hat.t$out_censor)
  # ===================================================================================
  # Gn.A1.t
  # ===================================================================================
  # plot initial fit
  if (verbose) lines(colMeans(Gn.A1.t) ~ T.uniq, col = 'yellow', lty = 1)
  # ===================================================================================
  # Qn.A1.t
  # ===================================================================================
  Qn.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

  # compute cumulative hazard
  # cum-product approach (2016-10-05)
  Qn.A1.t_full <- matrix(NA, nrow = n.data, ncol = ncol(h.hat.t_full))
  for (it in 1:n.data) {
    Qn.A1.t_full[it,] <- cumprod(1 - h.hat.t_full[it,])
  }
  Qn.A1.t <- Qn.A1.t_full[,T.uniq]

  # plot initial fit
  if (verbose) lines(colMeans(Qn.A1.t) ~ T.uniq, ...)
  # ===================================================================================
  # qn.A1.t
  # ===================================================================================
  # WILSON: rewrite in sweep?
  qn.A1.t_full <- matrix(0, nrow = n.data, ncol = ncol(Qn.A1.t_full))
  for (it.n in 1:n.data) {
    qn.A1.t_full[it.n,] <- h.hat.t_full[it.n,] * Qn.A1.t_full[it.n,]
  }
  qn.A1.t <- qn.A1.t_full[,T.uniq]

  # ===================================================================================
  # D1.t: calculate IC
  # D1.A1.t: calculate IC under intervention
  # ===================================================================================
  compute_IC <- function(dat, dW, T.uniq, h.hat.t_full, g.fitted, Gn.A1.t_full, Qn.A1.t, Qn.A1.t_full) {
    I.A.dW <- dat$A == dW
    n.data <- nrow(dat)

    D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
    D1.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

    for (it.n in 1:n.data) {

      t_Delta1.vec <- create_Yt_vector_with_censor(Time = dat$T.tilde[it.n], Delta = dat$delta[it.n], t.vec = 1:max(T.uniq))
      t.vec <- create_Yt_vector(Time = dat$T.tilde[it.n], t.vec = 1:max(T.uniq))
      alpha2 <- (t_Delta1.vec - t.vec * h.hat.t_full[it.n,])

      alpha1 <- -I.A.dW[it.n]/g.fitted[it.n]/Gn.A1.t_full[it.n,]/Qn.A1.t_full[it.n,]
      alpha1_A1 <- -1/g.fitted[it.n]/Gn.A1.t_full[it.n,]/Qn.A1.t_full[it.n,]

      not_complete <- alpha1 * alpha2
      not_complete_A1 <- alpha1_A1 * alpha2
      # D1 matrix
      D1.t[it.n, ] <- cumsum(not_complete)[T.uniq] * Qn.A1.t[it.n,] # complete influence curve
      D1.A1.t[it.n, ] <- cumsum(not_complete_A1)[T.uniq] * Qn.A1.t[it.n,] # also update those A = 0.
    }

    # turn unstable results to 0
    D1.t[is.na(D1.t)] <- 0
    D1.A1.t[is.na(D1.A1.t)] <- 0

    return(list(D1.t = D1.t,
                D1.A1.t = D1.A1.t))
  }

  initial_IC <- compute_IC(dat = dat,
                           dW = dW,
                           T.uniq = T.uniq,
                           h.hat.t_full = h.hat.t_full,
                           g.fitted = g.fitted,
                           Gn.A1.t_full = Gn.A1.t_full,
                           Qn.A1.t = Qn.A1.t,
                           Qn.A1.t_full = Qn.A1.t_full)
  D1.t <- initial_IC$D1.t
  D1.A1.t <- initial_IC$D1.A1.t
  # ===================================================================================
  # Pn.D1: efficient IC average
  # ===================================================================================
  # Pn.D1 vector
  Pn.D1.t <- colMeans(D1.t)
  # ===================================================================================
  # update
  # ===================================================================================
  message('targeting')
  stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq) # 10-17
  if(verbose) print(stopping.criteria)

  update.tensor <- matrix(0, nrow = n.data, ncol = length(T.uniq))
  iter.count <- 0
  stopping.prev <- Inf
  all_stopping <- numeric(stopping.criteria)
  all_loglikeli <- numeric()

  while ((stopping.criteria >= tol) & (iter.count <= max.iter)) { # ORGINAL
    # while ((stopping.criteria >= tol) & (iter.count <= max.iter) & ((stopping.prev - stopping.criteria) >= max(-tol, -1e-5))) { #WILSON: TEMPORARY
    if(verbose) print(stopping.criteria)
    # =============================================================================
    # update the qn
    # vectorized
    update.mat <- compute_onestep_update_matrix(D1.t.func.prev = D1.t,
                                                Pn.D1.func.prev = Pn.D1.t,
                                                dat = dat,
                                                T.uniq = T.uniq,
                                                W_names = W_names,
                                                dW = dW)
    # ------------------------------------------------------------------------
    update.tensor <- update.tensor + update.mat

    # accelerate when log-like becomes flat
    # if((stopping.prev - stopping.criteria) > 0 & (stopping.prev - stopping.criteria) < 1e-3) update.tensor <- update.tensor + update.mat*10

    # intergrand <- rowSums(update.tensor)
    # intergrand <- apply(update.tensor, c(1,2), sum)
    intergrand <- update.tensor
    intergrand[is.na(intergrand)] <- 0
    qn.current <- qn.A1.t * exp(epsilon.step * intergrand)
    qn.current_full <- qn.A1.t_full * exp(epsilon.step * replicate(T.max, intergrand[,1])) #10-23

    # For density sum > 1: normalize the updated qn
    norm.factor <- compute_step_cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf)[,1] #09-06
    qn.current[norm.factor > 1,] <- qn.current[norm.factor > 1,] / norm.factor[norm.factor > 1] #09-06
    qn.current_full[norm.factor > 1,] <- qn.current_full[norm.factor > 1,] / norm.factor[norm.factor > 1] #10-23

    # 11-26
    # For density sum > 1: truncate the density outside sum = 1 to be zero
    # i.e. flat cdf beyond sum to 1
    # cdf_per_subj <- compute_step_cdf(pdf.mat = qn.current, t.vec = T.uniq, start = -Inf)
    # qn.current[cdf_per_subj > 1] <- 0
    # cdf_per_subj <- compute_step_cdf(pdf.mat = qn.current_full, t.vec = 1:max(T.uniq), start = -Inf)
    # qn.current_full[cdf_per_subj > 1] <- 0

    # if some qn becomes all zero, prevent NA exisitence
    qn.current[is.na(qn.current)] <- 0
    qn.current_full[is.na(qn.current_full)] <- 0 #10-23
    # =============================================================================
    # compute new Qn
    Qn.current <- compute_step_cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf) # 2016-09-06
    cdf_offset <- 1 - Qn.current[,1] # 2016-09-06
    Qn.current <- Qn.current + cdf_offset # 2016-09-06

    Qn.current_full <- compute_step_cdf(pdf.mat = qn.current_full, t.vec = 1:max(T.uniq), start = Inf) # 10-23
    cdf_offset <- 1 - Qn.current_full[,1] # 10-23
    Qn.current_full <- Qn.current_full + cdf_offset # 10-23

    # check error
    # all.equal(compute_step_cdf(pdf.vec = qn.current[1,], t.vec = T.uniq, start = Inf), Qn.current[1,])
    # =============================================================================
    # compute new h_t
    h.hat.t_full_current <- matrix(0, nrow = n.data, ncol = max(T.uniq))
    for (it.n in 1:n.data) {
      h.hat.t_full_current[it.n, ] <- qn.current_full[it.n, ] / Qn.current_full[it.n,]
    }
    # compute new D1
    updated_IC <- compute_IC(dat = dat,
                             dW = dW,
                             T.uniq = T.uniq,
                             h.hat.t_full = h.hat.t_full_current,
                             g.fitted = g.fitted,
                             Gn.A1.t_full = Gn.A1.t_full,
                             Qn.A1.t = Qn.current,
                             Qn.A1.t_full = Qn.current_full)

    D1.t <- updated_IC$D1.t
    D1.A1.t <- updated_IC$D1.A1.t

    # compute new Pn.D1
    Pn.D1.t <- colMeans(D1.t)
    # ===================================================================================
    # previous stopping criteria
    stopping.prev <- stopping.criteria
    # new stopping criteria
    # stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq)
    stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq)/max(T.uniq))
    iter.count <- iter.count + 1
    # ===================================================================================
    # evaluate log-likelihood

    # construct obj
    obj <- list()
    obj$qn.current_full <- qn.current_full
    obj$Qn.current_full <- Qn.current_full
    obj$h.hat.t_full_current <- h.hat.t_full_current
    obj$dat <- dat
    # eval loglikeli
    loglike_here <- eval_loglike(obj, dW)
    all_loglikeli <- c(all_loglikeli, loglike_here)
    all_stopping <- c(all_stopping, stopping.criteria)


    ########################################################################
    # FOR DEBUG ONLY
    # if (TRUE) {
    #   # ------------------------------------------------------------
    #   # q <- seq(0,10,.1)
    #   ## truesurvExp <- 1 - pexp(q, rate = 1)
    #   # truesurvExp <- 1 - pexp(q, rate = .5)
    #   # plot(round(q*100,0), truesurvExp, type="l", cex=0.2, col = 'red', main = paste('l2 error =', stopping.criteria))
    #
    #   # library(survival)
    #   # n.data <- nrow(dat)
    #   # km.fit <- survfit(Surv(T,rep(1, n.data)) ~ A, data = dat)
    #   # lines(km.fit)
    #   # ------------------------------------------------------------
    #   # Psi.hat <- colMeans(Qn.current[dat$A==dW,])
    #
    #   # Psi.hat <- colMeans(Qn.current[dat$A==dW & dat$W==0,])
    #   # ------------------------------------------------------------
    #   # Q_weighted <- Qn.current/g.fitted[,1] # 09-18: inverse weight by propensity score
    #   # Q_weighted[dat$A!=dW,] <- 0
    #   # Psi.hat <- colMeans(Q_weighted) # 09-18: inverse weight by propensity score
    #   # ------------------------------------------------------------
    #     # 10-06: update all subjects with same W strata
    #   Psi.hat <- colMeans(Qn.current)
    #   # ------------------------------------------------------------
    #
    #   lines(Psi.hat ~ T.uniq, type = 'l', col = 'blue', lwd = .1)
    #   # ------------------------------------------------------------
    #   # legend('topright', lty=1, legend = c('true', 'KM', 'one-step'), col=c('red', 'black', 'blue'))
    # }
    ########################################################################
    if (iter.count == max.iter) {
      warning('Max Iter count reached, stop iteration.')
    }
  }

  if (!exists('Qn.current')) {
    # if the iteration immediately converge
    message('converge suddenly!')
    Qn.current <- Qn.A1.t
    updated_IC <- initial_IC
    Psi.hat <- colMeans(Qn.current)
  }

  # ===================================================================================
  # compute the target parameter
  # ===================================================================================
  # return the mean of those with observed A == dW
  Psi.hat <- colMeans(Qn.current)
  # variance of the EIC
  var_CI <- apply(updated_IC$D1.t, 2, var)/n.data
  # sup norm for each dim of EIC
  sup_norm_EIC <- abs(Pn.D1.t)

  variables <- list(T.uniq = T.uniq,
                    Qn.current = Qn.current,
                    D1.A1.t = D1.A1.t,
                    D1.t = D1.t,
                    Pn.D1.t = Pn.D1.t,
                    sup_norm_EIC = sup_norm_EIC)
  params <- list(stopping.criteria = stopping.criteria,
                 epsilon.step = epsilon.step,
                 iter.count = iter.count,
                 max.iter = max.iter,
                 dat = dat,
                 dW = dW)
  initial_fit <- list(h.hat.t = h.hat.t,
                      Qn.A1.t = Qn.A1.t,
                      qn.A1.t = qn.A1.t,
                      G.hat.t = G.hat.t,
                      g.fitted = g.fitted)
  convergence <- list(all_loglikeli = all_loglikeli,
                      all_stopping = all_stopping)
  # --------------------------------------------------
  to.return <- list(Psi.hat = Psi.hat,
                    T.uniq = T.uniq,
                    var = var_CI,
                    params = params,
                    variables = variables,
                    initial_fit = initial_fit,
                    convergence = convergence)
  class(to.return) <- 'surv_onestep'
  return(to.return)
}



#' One-step TMLE estimator for survival curve (No censoring)
#'
#' options to ADD:
#' SL.formula: the covariates to include in SL
#'
#' @param dat data.frame with columns T, A, W. All columns with character "W" will be treated as baseline covariates.
#' @param dW binary input vector specifying dynamic treatment (as a function output of W)
#' @param g.SL.Lib SuperLearner library for fitting treatment regression
#' @param ht.SL.Lib SuperLearner library for fitting conditional hazard regression
#' @param ... additional options for plotting initial fit curve
#' @param epsilon.step step size for one-step recursion
#' @param max.iter maximal number of recursion for one-step
#' @param T.cutoff  manual right censor the data; remove parts dont want to esimate
#' @param tol tolerance for optimization
#' @param verbose to plot the initial fit curve and the objective function value during optimzation
#'
#' @return Psi.hat vector of survival curve under intervention
#' @return T.uniq vector of time points where Psi.hat gets values (have same length as Psi.hat)
#' @return params list of meta-information of estimation
#' @return variables list of data summary
#' @return initial_fit list of initial fit (hazard, g_1, Delta)
#' @export
#'
#' @examples
#' @import dplyr
#' @import survtmle
#' @import abind
#' @import SuperLearner
surv_onestep_complete <- function(dat,
                                  dW,
                                  g.SL.Lib = c("SL.glm", "SL.step", "SL.glm.interaction"),
                                  ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                                  ...,
                                  epsilon.step = 1e-5, # the step size for one-step recursion
                                  max.iter = 1e3, # maximal number of recursion for one-step
                                  tol = 1/nrow(dat), # tolerance for optimization
                                  T.cutoff = NULL,
                                  verbose = TRUE) {
  # ================================================================================================
  # preparation
  # ================================================================================================
  after_check <- check_and_preprocess_data(dat = dat, dW = dW, T.cutoff = T.cutoff)
  dat <- after_check$dat
  dW <- after_check$dW
  n.data <- after_check$n.data
  W_names <- after_check$W_names

  W <- dat[,W_names]
  W <- as.data.frame(W)

  if(all(dW == 0)) {
    dat$A <- 1 - dat$A # when dW is all zero
    dW <- 1 - dW
  }else if(all(dW == 1)){

  }else{
    stop('not implemented!')
  }
  # ================================================================================================
  # estimate g(A|W)
  # ================================================================================================
  gHatSL <- SuperLearner(Y=dat$A, X=W, SL.library=g.SL.Lib, family="binomial")
  # g.hat for each observation
  g.fitted <- gHatSL$SL.predict
  # ================================================================================================
  # conditional hazard (by SL)
  # ================================================================================================
  message('estimating conditional hazard')
  T.uniq <- sort(unique(dat$T.tilde))
  T.max <- max(T.uniq)

  h.hat.t <- estimate_hazard_SL(dat = dat, T.uniq = T.uniq, ht.SL.Lib = ht.SL.Lib)
  # h.hat at all time t=[0,t.max]
  h.hat.t_full <- as.matrix(h.hat.t$out_haz_full)
  # h.hat at observed unique time t = T.grid
  h.hat.t <- as.matrix(h.hat.t$out_haz)
  # ================================================================================================
  # Qn.A1.t
  # ================================================================================================
  Qn.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

  # compute cumulative hazard
  # cum-product approach (2016-10-05)
  Qn.A1.t_full <- matrix(NA, nrow = n.data, ncol = ncol(h.hat.t_full))
  for (it in 1:n.data) {
    Qn.A1.t_full[it,] <- cumprod(1 - h.hat.t_full[it,])
  }
  Qn.A1.t <- Qn.A1.t_full[,T.uniq]

  # plot initial fit
  if (verbose) lines(colMeans(Qn.A1.t) ~ T.uniq, ...)
  # ================================================================================================
  # qn.A1.t
  # ================================================================================================
  # WILSON: rewrite in sweep?
  qn.A1.t_full <- matrix(0, nrow = n.data, ncol = ncol(Qn.A1.t_full))
  for (it.n in 1:n.data) {
    qn.A1.t_full[it.n,] <- h.hat.t_full[it.n,] * Qn.A1.t_full[it.n,]
  }
  qn.A1.t <- qn.A1.t_full[,T.uniq]

  # ================================================================================================
  # D1.t: calculate IC
  # D1.A1.t: calculate IC under intervention
  # ================================================================================================
  I.A.dW <- dat$A == dW

  D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
  D1.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

  for (it.n in 1:n.data) {
    Y.vec <- create_Yt_vector(Time = dat$T.tilde[it.n], t.vec = T.uniq)
    temp <- Y.vec - Qn.A1.t[it.n,]
    D1 <- temp / g.fitted[it.n] * I.A.dW[it.n]
    D1.A1 <- temp / g.fitted[it.n] # also update the samples without A = 1
    # D1 matrix
    D1.t[it.n,] <- D1
    D1.A1.t[it.n,] <- D1.A1
  }

  # ================================================================================================
  # Pn.D1: efficient IC average
  # ================================================================================================
  # Pn.D1 vector
  Pn.D1.t <- colMeans(D1.t)
  # ================================================================================================
  # update
  # ================================================================================================
  message('targeting')
  stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq) # 10-17

  update.tensor <- matrix(0, nrow = n.data, ncol = length(T.uniq))
  iter.count <- 0
  stopping.prev <- Inf

  # while ((stopping.criteria >= tol) & (iter.count <= max.iter)) { # ORGINAL
  while ((stopping.criteria >= tol) & (iter.count <= max.iter) & ((stopping.prev - stopping.criteria) >= max(-tol, -1e-5))) { #WILSON: TEMPORARY
    if(verbose) print(stopping.criteria)
    # =============================================================================
    # update the qn
    # vectorized
    update.mat <- compute_onestep_update_matrix(D1.t.func.prev = D1.t,
                                                Pn.D1.func.prev = Pn.D1.t,
                                                dat = dat,
                                                T.uniq = T.uniq,
                                                W_names = W_names,
                                                dW = dW)
    update.tensor <- update.tensor + update.mat

    # intergrand <- rowSums(update.tensor)
    # intergrand <- apply(update.tensor, c(1,2), sum)
    intergrand <- update.tensor
    intergrand[is.na(intergrand)] <- 0
    qn.current <- qn.A1.t * exp(epsilon.step * intergrand)

    # normalize the updated qn
    norm.factor <- compute_step_cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf)[,1] #09-06
    qn.current[norm.factor > 1,] <- qn.current[norm.factor > 1,] / norm.factor[norm.factor > 1] #09-06

    # 11-26
    # For density sum > 1: truncate the density outside sum = 1 to be zero
    # i.e. flat cdf beyond sum to 1
    # cdf_per_subj <- compute_step_cdf(pdf.mat = qn.current, t.vec = T.uniq, start = -Inf)
    # qn.current[cdf_per_subj > 1] <- 0

    # if some qn becomes all zero, prevent NA exisitence
    qn.current[is.na(qn.current)] <- 0
    # =============================================================================
    # compute new Qn

    Qn.current <- compute_step_cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf) # 2016-09-06
    cdf_offset <- 1 - Qn.current[,1] # 2016-09-06
    Qn.current <- Qn.current + cdf_offset # 2016-09-06

    # Qn.current <- apply(qn.current, 1, function(x) compute_step_cdf(pdf.vec = x, t.vec = T.uniq, start = Inf))
    # Qn.current <- t(Qn.current)

    # check error
    # all.equal(compute_step_cdf(pdf.vec = qn.current[1,], t.vec = T.uniq, start = Inf), Qn.current[1,])

    # < 2016-09-06
    # Qn.current <- matrix(NA, nrow = n.data, ncol = length(T.uniq))
    # for (it.n in 1:n.data) {
    # Qn.current[it.n,] <- rev(cumsum( rev(qn.current[it.n,]) ))
    # }

    # compute new D1
    D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
    D1.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
    for (it.n in 1:n.data) {
      Y.vec <- create_Yt_vector(Time = dat$T.tilde[it.n], t.vec = T.uniq)
      temp <- Y.vec - Qn.current[it.n,]
      D1 <- temp / g.fitted[it.n] * I.A.dW[it.n]
      D1.A1 <- temp / g.fitted[it.n] # also update the samples without A = 1
      # D1 matrix
      D1.t[it.n,] <- D1
      D1.A1.t[it.n,] <- D1.A1
    }
    # compute new Pn.D1
    Pn.D1.t <- colMeans(D1.t)
    # ================================================================================================
    # previous stopping criteria
    stopping.prev <- stopping.criteria
    # new stopping criteria
    stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq)
    iter.count <- iter.count + 1

    ########################################################################
    # FOR DEBUG ONLY
    # if (TRUE) {
    #   # ------------------------------------------------------------
    #   # q <- seq(0,10,.1)
    #   ## truesurvExp <- 1 - pexp(q, rate = 1)
    #   # truesurvExp <- 1 - pexp(q, rate = .5)
    #   # plot(round(q*100,0), truesurvExp, type="l", cex=0.2, col = 'red', main = paste('l2 error =', stopping.criteria))
    #
    #   # library(survival)
    #   # n.data <- nrow(dat)
    #   # km.fit <- survfit(Surv(T,rep(1, n.data)) ~ A, data = dat)
    #   # lines(km.fit)
    #   # ------------------------------------------------------------
    #   # Psi.hat <- colMeans(Qn.current[dat$A==dW,])
    #
    #   # Psi.hat <- colMeans(Qn.current[dat$A==dW & dat$W==0,])
    #   # ------------------------------------------------------------
    #   # Q_weighted <- Qn.current/g.fitted[,1] # 09-18: inverse weight by propensity score
    #   # Q_weighted[dat$A!=dW,] <- 0
    #   # Psi.hat <- colMeans(Q_weighted) # 09-18: inverse weight by propensity score
    #   # ------------------------------------------------------------
    #     # 10-06: update all subjects with same W strata
    #   Psi.hat <- colMeans(Qn.current)
    #   # ------------------------------------------------------------
    #
    #   lines(Psi.hat ~ T.uniq, type = 'l', col = 'blue', lwd = .1)
    #   # ------------------------------------------------------------
    #   # legend('topright', lty=1, legend = c('true', 'KM', 'one-step'), col=c('red', 'black', 'blue'))
    # }
    ########################################################################
    if (iter.count == max.iter) {
      warning('Max Iter count reached, stop iteration.')
    }
  }

  if (!exists('Qn.current')) {
    # if the iteration immediately converge
    message('converge suddenly!')
    Qn.current <- Qn.A1.t
    Psi.hat <- colMeans(Qn.current)
  }
  # ================================================================================================
  # compute the target parameter
  # ================================================================================================
  # return the mean of those with observed A == dW
  Psi.hat <- colMeans(Qn.current)
  # --------------------------------------------------
  variables <- list(T.uniq = T.uniq)
  params <- list(stopping.criteria = stopping.criteria,
                 epsilon.step = epsilon.step,
                 iter.count = iter.count,
                 max.iter = max.iter,
                 dat = dat)
  initial_fit <- list(h.hat.t = h.hat.t,
                      Qn.A1.t = Qn.A1.t,
                      qn.A1.t = qn.A1.t)
  to.return <- list(Psi.hat = Psi.hat,
                    T.uniq = T.uniq,
                    params = params,
                    variables = variables,
                    initial_fit = initial_fit)
  class(to.return) <- 'surv_onestep'
  return(to.return)
}


#' One-step TMLE estimator for survival curve
#'
#' one-step TMLE estimate of the difference of treatment-specific survival curves (S_1(t) - S_0(t)). Under right-censored data
#'
#' options to ADD:
#' SL.formula: the covariates to include in SL
#'
#' @param dat A data.frame with columns T.tilde, delta, A, W. T.tilde = min(T, C) is either the failure time of censor time, whichever happens first. 'delta'= I(T <= C) is the indicator of whether we observe failure time. A is binary treatment. W is baseline covariates. All columns with character "W" will be treated as baseline covariates.
#' @param dW A binary vector specifying dynamic treatment (as a function output of W)
#' @param g.SL.Lib A vector of string. SuperLearner library for fitting treatment regression
#' @param Delta.SL.Lib A vector of string. SuperLearner library for fitting censoring regression
#' @param ht.SL.Lib A vector of string. SuperLearner library for fitting conditional hazard regression
#' @param epsilon.step numeric. step size for one-step recursion
#' @param max.iter integer. maximal number of recursion for one-step
#' @param tol numeric. tolerance for optimization
#' @param T.cutoff int. Enforce randomized right-censoring to the observed data, so that don't estimate survival curve beyond a time point. Useful when time horizon is long.
#' @param verbose boolean. When TRUE, plot the initial fit curve, and output the objective function value during optimzation
#' @param ... additional options for plotting initial fit curve
#'
#' @return Psi.hat A numeric vector of estimated treatment-specific survival curve
#' @return T.uniq A vector of descrete time points where Psi.hat take values (have same length as Psi.hat)
#' @return params A list of estimation parameters set by user
#' @return variables A list of data summary statistics
#' @return initial_fit A list of initial fit values (hazard, g_1, Delta)
#'
#' @export
#'
#' @examples
#' library(simcausal)
#' D <- DAG.empty()
#'
#' D <- D +
#'     node("W", distr = "rbinom", size = 1, prob = .5) +
#'     node("A", distr = "rbinom", size = 1, prob = .15 + .5*W) +
#'     node("Trexp", distr = "rexp", rate = 1 + .5*W - .5*A) +
#'     node("Cweib", distr = "rweibull", shape = .7 - .2*W, scale = 1) +
#'     node("T", distr = "rconst", const = round(Trexp*100,0)) +
#'     node("C", distr = "rconst", const = round(Cweib*100, 0)) +
#'     node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
#'     node("delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
#' setD <- set.DAG(D)
#'
#' dat <- sim(setD, n=3e2)
#'
#' library(dplyr)
#' # only grab ID, W's, A, T.tilde, Delta
#' Wname <- grep('W', colnames(dat), value = TRUE)
#' dat <- dat[,c('ID', Wname, 'A', "T.tilde", "delta")]
#'
#' dW <- rep(1, nrow(dat))
#' onestepfit <- surv_onestep_difference(dat = dat,
#'                                        dW = dW,
#'                                        verbose = FALSE,
#'                                        epsilon.step = 1e-3,
#'                                        max.iter = 1e3)
#' @import dplyr
#' @import survtmle
#' @import abind
#' @import SuperLearner
surv_onestep_difference <- function(dat,
                                    dW = rep(1, nrow(dat)),
                                    g.SL.Lib = c("SL.glm", "SL.step", "SL.glm.interaction"),
                                    Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                                    ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                                    epsilon.step = 1e-5,
                                    max.iter = 1e3,
                                    tol = 1/nrow(dat),
                                    T.cutoff = NULL,
                                    verbose = TRUE,
                                    ...) {
  # ===================================================================================
  # preparation
  # ===================================================================================
  after_check <- check_and_preprocess_data(dat = dat, dW = dW, T.cutoff = T.cutoff)
  dat <- after_check$dat
  dW <- after_check$dW
  n.data <- after_check$n.data
  W_names <- after_check$W_names

  W <- dat[,W_names]
  W <- as.data.frame(W)

  # dW check
  dW = rep(1, nrow(dat))
  dat0 <- dat
  dat0$A <- 1 - dat0$A
  # if(all(dW == 0)) {
  #     dat$A <- 1 - dat$A # when dW is all zero
  #     dW <- 1 - dW
  # }else if(all(dW == 1)){

  # }else{
  #     stop('not implemented!')
  # }

  T.uniq <- sort(unique(dat$T.tilde))
  T.max <- max(T.uniq)
  # ===================================================================================
  # estimate g(A|W)
  # ===================================================================================
  gHatSL_1 <- SuperLearner(Y=dat$A, X=W, SL.library=g.SL.Lib, family="binomial")
  gHatSL_0 <- SuperLearner(Y=dat0$A, X=W, SL.library=g.SL.Lib, family="binomial")
  # g.hat for each observation
  g.fitted_1 <- gHatSL_1$SL.predict
  g.fitted_0 <- gHatSL_0$SL.predict
  # ===================================================================================
  # conditional hazard (by SL)
  # ===================================================================================
  message('estimating conditional hazard')

  h.hat.t_1 <- estimate_hazard_SL(dat = dat, T.uniq = T.uniq, ht.SL.Lib = ht.SL.Lib)
  h.hat.t_0 <- estimate_hazard_SL(dat = dat0, T.uniq = T.uniq, ht.SL.Lib = ht.SL.Lib)
  # h.hat at all time t=[0,t.max]
  h.hat.t_full_1 <- as.matrix(h.hat.t_1$out_haz_full)
  h.hat.t_full_0 <- as.matrix(h.hat.t_0$out_haz_full)
  # h.hat at observed unique time t = T.grid
  h.hat.t_1 <- as.matrix(h.hat.t_1$out_haz)
  h.hat.t_0 <- as.matrix(h.hat.t_0$out_haz)
  # ===================================================================================
  # estimate censoring G(A|W)
  # ===================================================================================
  message('estimating censoring')
  G.hat.t_1 <- estimate_censoring_SL(dat = dat, T.uniq = T.uniq,
                                     Delta.SL.Lib = Delta.SL.Lib)
  G.hat.t_0 <- estimate_censoring_SL(dat = dat0, T.uniq = T.uniq,
                                     Delta.SL.Lib = Delta.SL.Lib)
  # cutoff <- 0.1
  cutoff <- 0.05
  if(any(G.hat.t_1$out_censor_full <= cutoff)){
    warning('G.hat has extreme small values! lower truncate to 0.05')
    G.hat.t_1$out_censor_full[G.hat.t_1$out_censor_full < cutoff] <- cutoff
    G.hat.t_1$out_censor[G.hat.t_1$out_censor < cutoff] <- cutoff
  }
  if(any(G.hat.t_0$out_censor_full <= cutoff)){
    warning('G.hat has extreme small values! lower truncate to 0.05')
    G.hat.t_0$out_censor_full[G.hat.t_0$out_censor_full < cutoff] <- cutoff
    G.hat.t_0$out_censor[G.hat.t_0$out_censor < cutoff] <- cutoff
  }

  Gn.A1.t_full_1 <- as.matrix(G.hat.t_1$out_censor_full)
  Gn.A1.t_1 <- as.matrix(G.hat.t_1$out_censor)
  Gn.A1.t_full_0 <- as.matrix(G.hat.t_0$out_censor_full)
  Gn.A1.t_0 <- as.matrix(G.hat.t_0$out_censor)
  # ===================================================================================
  # Gn.A1.t
  # ===================================================================================
  # plot initial fit
  if (verbose) lines(colMeans(Gn.A1.t_1) ~ T.uniq, col = 'yellow', lty = 1)
  if (verbose) lines(colMeans(Gn.A1.t_0) ~ T.uniq, col = 'yellow', lty = 1)
  # ===================================================================================
  # Qn.A1.t
  # ===================================================================================
  Qn.A1.t_1 <- matrix(0, nrow = n.data, ncol = length(T.uniq))
  Qn.A1.t_0 <- matrix(0, nrow = n.data, ncol = length(T.uniq))

  # compute cumulative hazard
  # cum-product approach (2016-10-05)
  Qn.A1.t_full_1 <- matrix(NA, nrow = n.data, ncol = ncol(h.hat.t_full_1))
  Qn.A1.t_full_0 <- matrix(NA, nrow = n.data, ncol = ncol(h.hat.t_full_0))
  for (it in 1:n.data) {
    Qn.A1.t_full_1[it,] <- cumprod(1 - h.hat.t_full_1[it,])
  }
  Qn.A1.t_1 <- Qn.A1.t_full_1[,T.uniq]

  for (it in 1:n.data) {
    Qn.A1.t_full_0[it,] <- cumprod(1 - h.hat.t_full_0[it,])
  }
  Qn.A1.t_1 <- Qn.A1.t_full_1[,T.uniq]
  Qn.A1.t_0 <- Qn.A1.t_full_0[,T.uniq]

  # plot initial fit
  if (verbose) lines(colMeans(Qn.A1.t_1) ~ T.uniq)
  if (verbose) lines(colMeans(Qn.A1.t_0) ~ T.uniq)
  # ===================================================================================
  # qn.A1.t
  # ===================================================================================
  # WILSON: rewrite in sweep?
  qn.A1.t_full_1 <- matrix(0, nrow = n.data, ncol = ncol(Qn.A1.t_full_1))
  for (it.n in 1:n.data) {
    qn.A1.t_full_1[it.n,] <- h.hat.t_full_1[it.n,] * Qn.A1.t_full_1[it.n,]
  }
  qn.A1.t_1 <- qn.A1.t_full_1[,T.uniq]

  qn.A1.t_full_0 <- matrix(0, nrow = n.data, ncol = ncol(Qn.A1.t_full_0))
  for (it.n in 1:n.data) {
    qn.A1.t_full_0[it.n,] <- h.hat.t_full_0[it.n,] * Qn.A1.t_full_0[it.n,]
  }
  qn.A1.t_0 <- qn.A1.t_full_0[,T.uniq]

  # ===================================================================================
  # D1.t: calculate IC
  # D1.A1.t: calculate IC under intervention
  # ===================================================================================
  compute_IC <- function(dat, dW, T.uniq, h.hat.t_full, g.fitted, Gn.A1.t_full, Qn.A1.t, Qn.A1.t_full) {
    I.A.dW <- dat$A == dW
    n.data <- nrow(dat)

    D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
    D1.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

    for (it.n in 1:n.data) {

      t_Delta1.vec <- create_Yt_vector_with_censor(Time = dat$T.tilde[it.n], Delta = dat$delta[it.n], t.vec = 1:max(T.uniq))
      t.vec <- create_Yt_vector(Time = dat$T.tilde[it.n], t.vec = 1:max(T.uniq))
      alpha2 <- (t_Delta1.vec - t.vec * h.hat.t_full[it.n,])

      alpha1 <- -I.A.dW[it.n]/g.fitted[it.n]/Gn.A1.t_full[it.n,]/Qn.A1.t_full[it.n,]
      alpha1_A1 <- -1/g.fitted[it.n]/Gn.A1.t_full[it.n,]/Qn.A1.t_full[it.n,]

      not_complete <- alpha1 * alpha2
      not_complete_A1 <- alpha1_A1 * alpha2
      # D1 matrix
      D1.t[it.n, ] <- cumsum(not_complete)[T.uniq] * Qn.A1.t[it.n,] # complete influence curve
      D1.A1.t[it.n, ] <- cumsum(not_complete_A1)[T.uniq] * Qn.A1.t[it.n,] # also update those A = 0.
    }

    # turn unstable results to 0
    D1.t[is.na(D1.t)] <- 0
    D1.A1.t[is.na(D1.A1.t)] <- 0

    return(list(D1.t = D1.t,
                D1.A1.t = D1.A1.t))
  }

  initial_IC_1 <- compute_IC(dat = dat,
                             dW = rep(1, nrow(dat)),
                             T.uniq = T.uniq,
                             h.hat.t_full = h.hat.t_full_1,
                             g.fitted = g.fitted_1,
                             Gn.A1.t_full = Gn.A1.t_full_1,
                             Qn.A1.t = Qn.A1.t_1,
                             Qn.A1.t_full = Qn.A1.t_full_1)
  initial_IC_0 <- compute_IC(dat = dat0,
                             dW = rep(1, nrow(dat)),
                             T.uniq = T.uniq,
                             h.hat.t_full = h.hat.t_full_0,
                             g.fitted = g.fitted_0,
                             Gn.A1.t_full = Gn.A1.t_full_0,
                             Qn.A1.t = Qn.A1.t_0,
                             Qn.A1.t_full = Qn.A1.t_full_0)

  D1.t <- initial_IC_1$D1.t - initial_IC_0$D1.t
  D1.A1.t <- initial_IC_1$D1.A1.t - initial_IC_0$D1.A1.t
  # ===================================================================================
  # Pn.D1: efficient IC average
  # ===================================================================================
  # Pn.D1 vector
  Pn.D1.t <- colMeans(D1.t)
  # ===================================================================================
  # update
  # ===================================================================================
  message('targeting')
  stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq) # 10-17
  if(verbose) print(stopping.criteria)

  update.tensor <- matrix(0, nrow = n.data, ncol = length(T.uniq))
  iter.count <- 0
  stopping.prev <- Inf
  all_stopping <- numeric(stopping.criteria)
  all_loglikeli <- numeric()

  while ((stopping.criteria >= tol) & (iter.count <= max.iter)) { # ORGINAL
    # while ((stopping.criteria >= tol) & (iter.count <= max.iter) & ((stopping.prev - stopping.criteria) >= max(-tol, -1e-5))) { #WILSON: TEMPORARY
    if(verbose) print(stopping.criteria)
    # =============================================================================
    # update the qn
    # vectorized
    # update.mat <- compute_onestep_update_matrix_diff(D1.t.func.prev = D1.t,
    update.mat <- compute_onestep_update_matrix(D1.t.func.prev = D1.t,
                                                Pn.D1.func.prev = Pn.D1.t,
                                                dat = dat,
                                                T.uniq = T.uniq,
                                                W_names = W_names,
                                                dW = dW)
    update.tensor <- update.tensor + update.mat

    # accelerate when log-like becomes flat
    # if((stopping.prev - stopping.criteria) > 0 & (stopping.prev - stopping.criteria) < 1e-3) update.tensor <- update.tensor + update.mat*10

    # intergrand <- rowSums(update.tensor)
    # intergrand <- apply(update.tensor, c(1,2), sum)
    intergrand <- update.tensor
    intergrand[is.na(intergrand)] <- 0
    # qn.current_0 <- qn.A1.t_0 * exp(epsilon.step * intergrand)
    # qn.current_full_0 <- qn.A1.t_full_0 * exp(epsilon.step * replicate(T.max, intergrand[,1])) #10-23
    # qn.current_0 <- qn.A1.t_0 * exp(-epsilon.step * intergrand)
    # qn.current_full_0 <- qn.A1.t_full_0 * exp(-epsilon.step * replicate(T.max, intergrand[,1])) #10-23
    qn.current_0 <- qn.A1.t_0
    qn.current_full_0 <- qn.A1.t_full_0
    qn.current_1 <- qn.A1.t_1 * exp(epsilon.step * intergrand)
    qn.current_full_1 <- qn.A1.t_full_1 * exp(epsilon.step * replicate(T.max, intergrand[,1])) #10-23

    # For density sum > 1: normalize the updated qn
    norm.factor_1 <- compute_step_cdf(pdf.mat = qn.current_1, t.vec = T.uniq, start = Inf)[,1] #09-06
    # qn.current_1[norm.factor_1 > 1,] <- qn.current_1[norm.factor_1 > 1,] / norm.factor_1[norm.factor_1 > 1] #09-06
    # qn.current_full_1[norm.factor_1 > 1,] <- qn.current_full_1[norm.factor_1 > 1,] / norm.factor_1[norm.factor_1 > 1] #10-23
    norm.factor_0 <- compute_step_cdf(pdf.mat = qn.current_0, t.vec = T.uniq, start = Inf)[,1] #09-06
    # qn.current_0[norm.factor_0 > 1,] <- qn.current_0[norm.factor_0 > 1,] / norm.factor_0[norm.factor_0 > 1] #09-06
    # qn.current_full_0[norm.factor_0 > 1,] <- qn.current_full_0[norm.factor_0 > 1,] / norm.factor_0[norm.factor_0 > 1] #10-23

    # 11-26
    # For density sum > 1: truncate the density outside sum = 1 to be zero
    # i.e. flat cdf beyond sum to 1
    # cdf_per_subj <- compute_step_cdf(pdf.mat = qn.current, t.vec = T.uniq, start = -Inf)
    # qn.current[cdf_per_subj > 1] <- 0
    # cdf_per_subj <- compute_step_cdf(pdf.mat = qn.current_full, t.vec = 1:max(T.uniq), start = -Inf)
    # qn.current_full[cdf_per_subj > 1] <- 0

    # if some qn becomes all zero, prevent NA exisitence
    qn.current_0[is.na(qn.current_0)] <- 0
    qn.current_full_0[is.na(qn.current_full_0)] <- 0 #10-23
    qn.current_1[is.na(qn.current_1)] <- 0
    qn.current_full_1[is.na(qn.current_full_1)] <- 0 #10-23
    # =============================================================================
    # compute new Qn

    Qn.current_1 <- compute_step_cdf(pdf.mat = qn.current_1, t.vec = T.uniq, start = Inf) # 2016-09-06
    cdf_offset_1 <- 1 - Qn.current_1[,1] # 2016-09-06
    Qn.current_1 <- Qn.current_1 + cdf_offset_1 # 2016-09-06

    Qn.current_full_1 <- compute_step_cdf(pdf.mat = qn.current_full_1, t.vec = 1:max(T.uniq), start = Inf) # 10-23
    cdf_offset_1 <- 1 - Qn.current_full_1[,1] # 10-23
    Qn.current_full_1 <- Qn.current_full_1 + cdf_offset_1 # 10-23

    Qn.current_0 <- compute_step_cdf(pdf.mat = qn.current_0, t.vec = T.uniq, start = Inf) # 2016-09-06
    cdf_offset_0 <- 1 - Qn.current_0[,1] # 2016-09-06
    Qn.current_0 <- Qn.current_0 + cdf_offset_0 # 2016-09-06

    Qn.current_full_0 <- compute_step_cdf(pdf.mat = qn.current_full_0, t.vec = 1:max(T.uniq), start = Inf) # 10-23
    cdf_offset_0 <- 1 - Qn.current_full_0[,1] # 10-23
    Qn.current_full_0 <- Qn.current_full_0 + cdf_offset_0 # 10-23

    Psin.current <- Qn.current_1 - Qn.current_0
    # check error
    # all.equal(compute_step_cdf(pdf.vec = qn.current[1,], t.vec = T.uniq, start = Inf), Qn.current[1,])
    # =============================================================================
    # compute new h_t
    h.hat.t_full_current_1 <- matrix(0, nrow = n.data, ncol = max(T.uniq))
    h.hat.t_full_current_0 <- matrix(0, nrow = n.data, ncol = max(T.uniq))
    for (it.n in 1:n.data) {
      h.hat.t_full_current_1[it.n, ] <- qn.current_full_1[it.n, ] / Qn.current_full_1[it.n,]
    }
    for (it.n in 1:n.data) {
      h.hat.t_full_current_0[it.n, ] <- qn.current_full_0[it.n, ] / Qn.current_full_0[it.n,]
    }
    # compute new D1
    updated_IC_1 <- compute_IC(dat = dat,
                               dW = 1,
                               T.uniq = T.uniq,
                               h.hat.t_full = h.hat.t_full_current_1,
                               g.fitted = g.fitted_1,
                               Gn.A1.t_full = Gn.A1.t_full_1,
                               Qn.A1.t = Qn.current_1,
                               Qn.A1.t_full = Qn.current_full_1)


    updated_IC_0 <- compute_IC(dat = dat0,
                               dW = 1,
                               T.uniq = T.uniq,
                               h.hat.t_full = h.hat.t_full_current_0,
                               g.fitted = g.fitted_0,
                               Gn.A1.t_full = Gn.A1.t_full_0,
                               Qn.A1.t = Qn.current_0,
                               Qn.A1.t_full = Qn.current_full_0)

    D1.t <- updated_IC_1$D1.t - updated_IC_0$D1.t
    D1.A1.t <- updated_IC_1$D1.A1.t - updated_IC_0$D1.A1.t

    # compute new Pn.D1
    Pn.D1.t <- colMeans(D1.t)
    # ===================================================================================
    # previous stopping criteria
    stopping.prev <- stopping.criteria
    # new stopping criteria
    # stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq)
    stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq)/max(T.uniq))
    iter.count <- iter.count + 1
    # ===================================================================================
    # evaluate log-likelihood

    # construct obj
    # obj <- list()
    # obj$qn.current_full_1 <- qn.current_full_1
    # obj$Qn.current_full_1 <- Qn.current_full_1
    # obj$h.hat.t_full_current_1 <- h.hat.t_full_current_1
    # obj$dat <- dat
    # obj$qn.current_full_0 <- qn.current_full_0
    # obj$Qn.current_full_0 <- Qn.current_full_0
    # obj$h.hat.t_full_current_0 <- h.hat.t_full_current_0
    # obj$dat <- dat

    # obj$Psin.current <- Psin.current

    # eval loglikeli
    # loglike_here <- eval_loglike(obj, dW)
    # all_loglikeli <- c(all_loglikeli, loglike_here)
    # all_stopping <- c(all_stopping, stopping.criteria)


    ########################################################################
    # FOR DEBUG ONLY
    # if (TRUE) {
    #   # ------------------------------------------------------------
    #   # q <- seq(0,10,.1)
    #   ## truesurvExp <- 1 - pexp(q, rate = 1)
    #   # truesurvExp <- 1 - pexp(q, rate = .5)
    #   # plot(round(q*100,0), truesurvExp, type="l", cex=0.2, col = 'red', main = paste('l2 error =', stopping.criteria))
    #
    #   # library(survival)
    #   # n.data <- nrow(dat)
    #   # km.fit <- survfit(Surv(T,rep(1, n.data)) ~ A, data = dat)
    #   # lines(km.fit)
    #   # ------------------------------------------------------------
    #   # Psi.hat <- colMeans(Qn.current[dat$A==dW,])
    #
    #   # Psi.hat <- colMeans(Qn.current[dat$A==dW & dat$W==0,])
    #   # ------------------------------------------------------------
    #   # Q_weighted <- Qn.current/g.fitted[,1] # 09-18: inverse weight by propensity score
    #   # Q_weighted[dat$A!=dW,] <- 0
    #   # Psi.hat <- colMeans(Q_weighted) # 09-18: inverse weight by propensity score
    #   # ------------------------------------------------------------
    #     # 10-06: update all subjects with same W strata
    #   # Psi.hat <- colMeans(Qn.current)
    #   # ------------------------------------------------------------
    #   lines(colMeans(Qn.current_1) ~ T.uniq)
    #   lines(colMeans(Qn.current_0) ~ T.uniq)
    #   # lines(Psi.hat ~ T.uniq, type = 'l', col = 'blue', lwd = .1)
    #   lines(colMeans(Psin.current) ~ T.uniq, type = 'l', col = 'blue', lwd = .1)
    #   # ------------------------------------------------------------
    #   # legend('topright', lty=1, legend = c('true', 'KM', 'one-step'), col=c('red', 'black', 'blue'))
    # }
    ########################################################################
    if (iter.count == max.iter) {
      warning('Max Iter count reached, stop iteration.')
    }
  }

  if (!exists('Qn.current_1')) {
    # if the iteration immediately converge
    message('converge suddenly!')
    Qn.current <- Qn.A1.t_1 - Qn.A1.t_0
    updated_IC <- initial_IC_1 - initial_IC_0
    Psi.hat <- colMeans(Qn.current)
  }

  # ===================================================================================
  # compute the target parameter
  # ===================================================================================
  # return the mean of those with observed A == dW
  # Psi.hat <- colMeans(Qn.current)
  Psi.hat <- colMeans(Psin.current)
  # variance of the EIC
  # var_CI <- apply(updated_IC$D1.t, 2, var)/n.data
  var_CI <- apply(D1.t, 2, var)/n.data
  # --------------------------------------------------
  # sup norm for each dim of EIC
  sup_norm_EIC <- abs(Pn.D1.t)

  variables <- list(T.uniq = T.uniq,
                    Psin.current = Psin.current,
                    D1.A1.t = D1.A1.t,
                    D1.t = D1.t,
                    Pn.D1.t = Pn.D1.t,
                    sup_norm_EIC = sup_norm_EIC)
  params <- list(stopping.criteria = stopping.criteria,
                 epsilon.step = epsilon.step,
                 iter.count = iter.count,
                 max.iter = max.iter,
                 dat = dat,
                 dW = dW)
  initial_fit_1 <- list(h.hat.t_1 = h.hat.t_1,
                        Qn.A1.t_1 = Qn.A1.t_1,
                        qn.A1.t_1 = qn.A1.t_1,
                        G.hat.t_1 = G.hat.t_1)
  initial_fit_0 <- list(h.hat.t_0 = h.hat.t_0,
                        Qn.A1.t_0 = Qn.A1.t_0,
                        qn.A1.t_0 = qn.A1.t_0,
                        G.hat.t_0 = G.hat.t_0)
  convergence <- list(all_loglikeli = all_loglikeli,
                      all_stopping = all_stopping)
  # --------------------------------------------------
  to.return <- list(Psi.hat = Psi.hat,
                    T.uniq = T.uniq,
                    var = var_CI,
                    params = params,
                    variables = variables,
                    initial_fit_1 = initial_fit_1,
                    initial_fit_0 = initial_fit_0,
                    convergence = convergence)
  class(to.return) <- 'surv_onestep'
  return(to.return)

}


#' Perform one-step TMLE update of survival curve
#'
#' @param D1.t.func.prev n*p matrix of previous influence curve
#' @param Pn.D1.func.prev p vector of previous mean influence curve
#' @param dat input data.frame
#' @param T.uniq grid of unique event times
#' @param W_names vector of the names of baseline covariates
#' @param dW dynamic intervention
#'
#' @return
#' @export
#'
#' @examples
#' # TO DO
#' @importFrom dplyr left_join
compute_onestep_update_matrix <- function(D1.t.func.prev, Pn.D1.func.prev, dat, T.uniq, W_names, dW) {
  # calculate the number inside exp{} expression in submodel
  # each strata of Q is updated the same
  browser()
  numerator1 <- l2_inner_prod_step(Pn.D1.func.prev, D1.t.func.prev, T.uniq)

  numerator <- sweep(D1.t.func.prev, MARGIN=2, -abs(Pn.D1.func.prev),`*`)
  # numerator <- sweep(D1.t.func.prev, MARGIN=2, abs(Pn.D1.func.prev),`*`) # WROOOOOONG
  result <- numerator /
    sqrt(l2_inner_prod_step(Pn.D1.func.prev, Pn.D1.func.prev, T.uniq))

  return(result)
}

# =======================================================================================
# one-step TMLE for survival at a specific end-point
# =======================================================================================

#' One-step TMLE estimator for survival at specific time point
#'
#' @param dat data.frame with columns T, A, C, W. All columns with character "W" will be treated as baseline covariates.
#' @param tk time point to compute survival probability
#' @param dW binary input vector specifying dynamic treatment (as a function output of W)
#' @param SL.trt SuperLearner library for fitting treatment regression
#' @param SL.ctime SuperLearner library for fitting censoring regression
#' @param SL.ftime SuperLearner library for fitting conditional hazard regression
#' @param maxIter maximal number of recursion for one-step
#' @param epsilon_step step size for one-step recursion
#' @param tol tolerance for optimization
#' @param T.cutoff  manual right censor the data; remove parts dont want to esimate
#' @param verbose to print log-likelihood value during optimzation
#'
#' @return
#' @export
#'
#' @examples
#' # TO DO
#' @import dplyr
#' @import survtmle
onestep_single_t <- function(dat, tk, dW = rep(1, nrow(dat)),
                             SL.trt = c("SL.glm", "SL.step", "SL.earth"),
                             SL.ctime = c("SL.glm", "SL.step", "SL.earth"),
                             SL.ftime = c("SL.glm", "SL.step", "SL.earth"),
                             maxIter = 3e2,
                             epsilon_step = 1e-3,
                             tol = 1/nrow(dat),
                             T.cutoff = NULL,
                             verbose = FALSE){
  # ====================================================================================================
  # input validation
  # ====================================================================================================
  after_check <- check_and_preprocess_data(dat = dat, dW = dW, T.cutoff = T.cutoff)
  dat <- after_check$dat
  dW <- after_check$dW
  n.data <- after_check$n.data
  W_names <- after_check$W_names
  # ====================================================================================================
  # preparation: make data in survtmle format (dat_david)
  # ====================================================================================================
  # transform original data into SL-friendly format
  dat_david <- dat

  dat_david <- rename(dat_david, ftime = T.tilde)
  dat_david <- rename(dat_david, trt = A)

  if(all(dW == 0)) {
    dat_david$trt <- 1 - dat_david$trt # when dW is all zero
  }else if(all(dW == 1)){

  }else{
    stop('not implemented!')
  }

  if ('ID' %in% toupper(colnames(dat_david))) {
    # if there are already id in the dataset
    dat_david <- rename(dat_david, id = ID)
  }else{
    warning('no id exist, create \'id\' on our own')
    dat_david$id <- 1:nrow(dat_david)
  }

  # censoring
  dat_david <- rename(dat_david, ftype = delta)

  # remove all other useless columns
  baseline_name <- W_names
  keeps <- c("id", baseline_name, 'ftime', 'ftype', 'trt')
  dat_david <- dat_david[,keeps]
  # ====================================================================================================
  # prepare
  # ====================================================================================================
  T.uniq <- unique(sort(dat_david$ftime))
  T.max <- max(T.uniq)

  adjustVars <- dat_david[,baseline_name, drop = FALSE]
  # ====================================================================================================
  # estimate g
  # ====================================================================================================
  message('estimating g_1')
  g1_hat <- estimateTreatment(dat = dat_david, adjustVars = adjustVars,
                              SL.trt = SL.trt,verbose = verbose, returnModels = FALSE)
  g1_dat <- g1_hat$dat
  # ====================================================================================================
  # make datalist
  # ====================================================================================================
  datalist <- survtmle::makeDataList(dat = g1_dat,
                                     J = 1, # one kind of failure
                                     ntrt = 2, # one kind of treatment
                                     uniqtrt = c(0,1),
                                     t0 = tk, # time to predict on
                                     bounds=NULL)
  # yo <- datalist[[3]]
  # ====================================================================================================
  # estimate g_2 (censoring)
  # ====================================================================================================
  message('estimating g_2')
  g2_hat <- survtmle:::estimateCensoring(dataList = datalist, adjustVars = adjustVars,
                                         t0 = tk,
                                         ntrt = 2, # one kind of treatment
                                         uniqtrt = c(0,1),
                                         SL.ctime = SL.ctime,
                                         returnModels = FALSE,verbose = verbose)
  dataList2 <- g2_hat$dataList
  # ====================================================================================================
  # estimate h(t) (hazard)
  # ====================================================================================================
  message('estimating hazard')
  h_hat <- survtmle::estimateHazards(dataList = dataList2,
                                     J = 1,
                                     adjustVars = adjustVars,
                                     SL.ftime = SL.ftime,
                                     returnModels = FALSE,
                                     verbose = verbose,
                                     glm.ftime = NULL)
  dataList2 <- h_hat$dataList
  # check convergence
  suppressWarnings(if (all(dataList2[[1]] == "convergence failure")) {
    return("estimation convergence failure")
  })
  # ====================================================================================================
  # transform to survivial
  # ====================================================================================================
  dataList2 <- updateVariables(dataList = dataList2, allJ = 1,
                               ofInterestJ = 1, nJ = 2, uniqtrt = c(0,1),
                               ntrt = 2, t0 = tk, verbose = verbose)
  # hehe <- dataList2[[3]]
  # ====================================================================================================
  # get IC
  # ====================================================================================================
  dat_david2 <- getHazardInfluenceCurve(dataList = dataList2, dat = dat_david,
                                        ofInterestJ = 1, allJ = 1, nJ = 2, uniqtrt = c(0,1),
                                        ntrt = 2, verbose = verbose, t0 = tk)
  infCurves <- dat_david2[, grep("D.j", names(dat_david2))]
  meanIC <- colMeans(infCurves)
  # ====================================================================================================
  # targeting
  # ====================================================================================================
  calcLoss <- function(Y, QAW){
    -mean(Y * log(QAW) + (1-Y) * log(1 - QAW))
  }


  # if the derivative of the loss the positive, change the tergeting direction
  epsilon_step1 <- epsilon_step2 <- epsilon_step
  if (meanIC[2] < 0) { epsilon_step2 <- -epsilon_step2}
  if (meanIC[1] < 0) { epsilon_step1 <- -epsilon_step1}

  loss_old <- Inf
  # loss_new <- calcLoss(Y = dataList2$`1`$N1, QAW = dataList2$`1`$Q1Haz)
  loss_new <- calcLoss(Y = dataList2$obs$N1, QAW = dataList2$obs$Q1Haz)
  message('targeting')
  iter_count <- 0

  # while (any(abs(meanIC) > tol) & iter_count <= maxIter) {
  # while (any(abs(meanIC[2]) > tol) & iter_count <= maxIter) {
  while ((loss_new <= loss_old) & iter_count <= maxIter) {
    iter_count <- iter_count + 1
    # print(loss_new)
    # print(meanIC[1,])

    # fluctuate -> update to dataList2
    dataList2$`1`$Q1Haz <- plogis(qlogis(dataList2$`1`$Q1Haz) + epsilon_step2 * dataList2$`1`$H1.jSelf.z1 + epsilon_step1 * dataList2$`1`$H1.jSelf.z0)
    dataList2$`0`$Q1Haz <- plogis(qlogis(dataList2$`0`$Q1Haz) + epsilon_step2 * dataList2$`0`$H1.jSelf.z1 + epsilon_step1 * dataList2$`0`$H1.jSelf.z0)
    dataList2$obs$Q1Haz <- plogis(qlogis(dataList2$obs$Q1Haz) + epsilon_step2 * dataList2$obs$H1.jSelf.z1 + epsilon_step1 * dataList2$obs$H1.jSelf.z0)


    # calculate survival again
    dataList2 <- updateVariables(dataList = dataList2, allJ = 1,
                                 ofInterestJ = 1, nJ = 2, uniqtrt = c(0,1),
                                 ntrt = 2, t0 = tk, verbose = verbose)
    # calculate IC again
    dat_david2 <- getHazardInfluenceCurve(dataList = dataList2, dat = dat_david2,
                                          ofInterestJ = 1, allJ = 1, nJ = 2, uniqtrt = c(0,1),
                                          ntrt = 2, verbose = verbose, t0 = tk)
    infCurves <- dat_david2[, grep("D.j", names(dat_david2))]
    meanIC_old <- meanIC
    meanIC <- colMeans(infCurves)

    # loss_new <- calcLoss(Y = dataList2$`1`$N1, QAW = dataList2$`1`$Q1Haz)
    loss_old <- loss_new
    loss_new <- calcLoss(Y = dataList2$obs$N1, QAW = dataList2$obs$Q1Haz)


    # if one converges, then stop update
    if ((abs(meanIC[1]) < tol) | (meanIC_old[1] * meanIC[1] <= 0)) {
      # if changes sign or converges, then stop update
      epsilon_step1 <- 0
    }
    if ((abs(meanIC[2]) < tol) | (meanIC_old[2] * meanIC[2] <= 0)) {
      # if changes sign or converges, then stop update
      epsilon_step2 <- 0
    }
    # if (all(abs(meanIC) < tol)) {
    if (abs(meanIC)[2] < tol) {
      # all IC become zero mean
      message('Success Converge!')
      break()
    }
    if (epsilon_step1 == 0 & epsilon_step2 == 0){
      message('Success Converge! with epsilon_step too large.')
      break()
    }
  }

  if (iter_count == maxIter + 1) {
    warning("TMLE fluctuations did not converge. Check that meanIC is adequately small and proceed with caution.")
  }
  # ====================================================================================================
  # get final estimates
  # ====================================================================================================
  est <- rowNames <- NULL

  # parameter estimates
  for (j in 1) {
    for (z in c(0,1)) {
      eval(parse(text = paste("est <- rbind(est, dat_david2$margF",
                              j, ".z", z, ".t0[1])", sep = "")))
      rowNames <- c(rowNames, paste(c(z, j), collapse = " "))
    }
  }
  row.names(est) <- rowNames
  var <- t(as.matrix(infCurves)) %*% as.matrix(infCurves)/n.data^2
  row.names(var) <- colnames(var) <- rowNames

  # output static interventions
  est <- 1 - est['1 1',]
  var <- var['1 1', '1 1']

  # if (all(dW == 1)) {
  #     est <- 1 - est['1 1',]
  #     var <- var['1 1', '1 1']
  # }else if (all(dW == 0)) {
  #     est <- 1 - est['0 1',]
  #     var <- var['0 1', '0 1']
  # }

  return(list(est = est, var = var, meanIC = meanIC, ic = infCurves))
}


#' One-step TMLE estimator for survival at specific time point; Loop over all times
#'
#' @param dat data.frame with columns T, A, C, W. All columns with character "W" will be treated as baseline covariates.
#' @param tk time point to compute survival probability
#' @param dW binary input vector specifying dynamic treatment (as a function output of W)
#' @param SL.trt SuperLearner library for fitting treatment regression
#' @param SL.ctime SuperLearner library for fitting censoring regression
#' @param SL.ftime SuperLearner library for fitting conditional hazard regression
#' @param maxIter maximal number of recursion for one-step
#' @param epsilon_step step size for one-step recursion
#' @param tol tolerance for optimization
#' @param T.cutoff  manual right censor the data; remove parts dont want to esimate
#' @param verbose to print log-likelihood value during optimzation
#'
#' @return
#' @export
#'
#' @examples
#' # TO DO
#' @import dplyr
#' @import survtmle
onestep_single_t_loopall <- function(dat, dW = rep(1, nrow(dat)),
                                     SL.trt = c("SL.glm", "SL.step", "SL.earth"),
                                     SL.ctime = c("SL.glm", "SL.step", "SL.earth"),
                                     SL.ftime = c("SL.glm", "SL.step", "SL.earth"),
                                     maxIter = 3e2,
                                     epsilon_step = 1e-3,
                                     tol = 1/nrow(dat),
                                     T.cutoff = NULL,
                                     verbose = FALSE){
  # ====================================================================================================
  # input validation
  # ====================================================================================================
  after_check <- check_and_preprocess_data(dat = dat, dW = dW, T.cutoff = T.cutoff)
  dat <- after_check$dat
  dW <- after_check$dW
  n.data <- after_check$n.data
  W_names <- after_check$W_names
  # ====================================================================================================
  # preparation: make data in survtmle format (dat_david)
  # ====================================================================================================
  # transform original data into SL-friendly format
  dat_david <- dat

  dat_david <- rename(dat_david, ftime = T.tilde)
  dat_david <- rename(dat_david, trt = A)

  if(all(dW == 0)) {
    dat_david$trt <- 1 - dat_david$trt # when dW is all zero
  }else if(all(dW == 1)){

  }else{
    stop('not implemented!')
  }

  if ('ID' %in% toupper(colnames(dat_david))) {
    # if there are already id in the dataset
    dat_david <- rename(dat_david, id = ID)
  }else{
    warning('no id exist, create \'id\' on our own')
    dat_david$id <- 1:nrow(dat_david)
  }

  # censoring
  dat_david <- rename(dat_david, ftype = delta)

  # remove all other useless columns
  baseline_name <- W_names
  keeps <- c("id", baseline_name, 'ftime', 'ftype', 'trt')
  dat_david <- dat_david[,keeps]
  # ====================================================================================================
  # prepare
  # ====================================================================================================
  T.uniq <- unique(sort(dat_david$ftime))
  T.max <- max(T.uniq)

  adjustVars <- dat_david[,baseline_name, drop = FALSE]
  # ====================================================================================================
  # estimate g
  # ====================================================================================================
  message('estimating g_1')
  g1_hat <- estimateTreatment(dat = dat_david, adjustVars = adjustVars,
                              SL.trt = SL.trt,verbose = verbose, returnModels = FALSE)
  g1_dat <- g1_hat$dat
  # ====================================================================================================
  # make datalist
  # ====================================================================================================
  datalist <- survtmle::makeDataList(dat = g1_dat,
                                     J = 1, # one kind of failure
                                     ntrt = 2, # one kind of treatment
                                     uniqtrt = c(0,1),
                                     t0 = T.max, # time to predict on
                                     bounds=NULL)
  # yo <- datalist[[3]]
  # ====================================================================================================
  # estimate g_2 (censoring)
  # ====================================================================================================
  message('estimating g_2')
  g2_hat <- survtmle:::estimateCensoring(dataList = datalist, adjustVars = adjustVars,
                                         t0 = T.max,
                                         ntrt = 2, # one kind of treatment
                                         uniqtrt = c(0,1),
                                         SL.ctime = SL.ctime,
                                         returnModels = FALSE,verbose = verbose)
  dataList2 <- g2_hat$dataList
  # ====================================================================================================
  # estimate h(t) (hazard)
  # ====================================================================================================
  message('estimating hazard')
  h_hat <- survtmle::estimateHazards(dataList = dataList2,
                                     J = 1,
                                     adjustVars = adjustVars,
                                     SL.ftime = SL.ftime,
                                     returnModels = FALSE,
                                     verbose = verbose,
                                     glm.ftime = NULL)
  dataList2 <- h_hat$dataList
  # check convergence
  suppressWarnings(if (all(dataList2[[1]] == "convergence failure")) {
    return("estimation convergence failure")
  })
  # ====================================================================================================
  # transform to survivial
  # ====================================================================================================
  dataList2 <- updateVariables(dataList = dataList2, allJ = 1,
                               ofInterestJ = 1, nJ = 2, uniqtrt = c(0,1),
                               ntrt = 2, t0 = T.max, verbose = verbose)
  # hehe <- dataList2[[3]]
  dataList2_before_target <- dataList2
  onestep_out_all <- list()
  tk_count <- 0
  for (tk in T.uniq) {
    tk_count <- tk_count + 1
    dataList2 <- dataList2_before_target

    # ====================================================================================================
    # get IC
    # ====================================================================================================
    dat_david2 <- getHazardInfluenceCurve(dataList = dataList2, dat = dat_david,
                                          ofInterestJ = 1, allJ = 1, nJ = 2, uniqtrt = c(0,1),
                                          ntrt = 2, verbose = verbose, t0 = tk)
    infCurves <- dat_david2[, grep("D.j", names(dat_david2))]
    meanIC <- colMeans(infCurves)
    # ====================================================================================================
    # targeting
    # ====================================================================================================
    calcLoss <- function(Y, QAW){
      -mean(Y * log(QAW) + (1-Y) * log(1 - QAW))
    }

    # if the derivative of the loss the positive, change the tergeting direction
    epsilon_step1 <- epsilon_step2 <- epsilon_step
    if (meanIC[2] < 0) { epsilon_step2 <- -epsilon_step2}
    if (meanIC[1] < 0) { epsilon_step1 <- -epsilon_step1}

    loss_old <- Inf
    # loss_new <- calcLoss(Y = dataList2$`1`$N1, QAW = dataList2$`1`$Q1Haz)
    loss_new <- calcLoss(Y = dataList2$obs$N1, QAW = dataList2$obs$Q1Haz)
    message(paste('targeting', tk))
    iter_count <- 0

    # while (any(abs(meanIC) > tol) & iter_count <= maxIter) {
    # while (any(abs(meanIC[2]) > tol) & iter_count <= maxIter) {
    while ((loss_new <= loss_old) & iter_count <= maxIter) {
      iter_count <- iter_count + 1
      if(verbose) print(loss_new)
      # print(meanIC[1,])

      # fluctuate -> update to dataList2
      dataList2$`1`$Q1Haz <- plogis(qlogis(dataList2$`1`$Q1Haz) + epsilon_step2 * dataList2$`1`$H1.jSelf.z1 + epsilon_step1 * dataList2$`1`$H1.jSelf.z0)
      dataList2$`0`$Q1Haz <- plogis(qlogis(dataList2$`0`$Q1Haz) + epsilon_step2 * dataList2$`0`$H1.jSelf.z1 + epsilon_step1 * dataList2$`0`$H1.jSelf.z0)
      dataList2$obs$Q1Haz <- plogis(qlogis(dataList2$obs$Q1Haz) + epsilon_step2 * dataList2$obs$H1.jSelf.z1 + epsilon_step1 * dataList2$obs$H1.jSelf.z0)


      # calculate survival again
      dataList2 <- updateVariables(dataList = dataList2, allJ = 1,
                                   ofInterestJ = 1, nJ = 2, uniqtrt = c(0,1),
                                   ntrt = 2, t0 = tk, verbose = verbose)
      # calculate IC again
      dat_david2 <- getHazardInfluenceCurve(dataList = dataList2, dat = dat_david2,
                                            ofInterestJ = 1, allJ = 1, nJ = 2, uniqtrt = c(0,1),
                                            ntrt = 2, verbose = verbose, t0 = tk)
      infCurves <- dat_david2[, grep("D.j", names(dat_david2))]
      meanIC_old <- meanIC
      meanIC <- colMeans(infCurves)

      # loss_new <- calcLoss(Y = dataList2$`1`$N1, QAW = dataList2$`1`$Q1Haz)
      loss_old <- loss_new
      loss_new <- calcLoss(Y = dataList2$obs$N1, QAW = dataList2$obs$Q1Haz)


      # if one converges, then stop update
      if ((abs(meanIC[1]) < tol) | (meanIC_old[1] * meanIC[1] <= 0)) {
        # if changes sign or converges, then stop update
        epsilon_step1 <- 0
      }
      if ((abs(meanIC[2]) < tol) | (meanIC_old[2] * meanIC[2] <= 0)) {
        # if changes sign or converges, then stop update
        epsilon_step2 <- 0
      }

      # if (all(abs(meanIC) < tol)) {
      # print(abs(meanIC))
      if (abs(meanIC)[2] < tol) {
        # all IC become zero mean
        message('Success Converge!')
        break()
      }
      if (epsilon_step1 == 0 & epsilon_step2 == 0){
        message('Success Converge! with epsilon_step too large.')
        break()
      }
    }

    if (iter_count == maxIter + 1) {
      message("TMLE fluctuations did not converge. Check that meanIC is adequately small and proceed with caution.")
    }


    # ====================================================================================================
    # get final estimates
    # ====================================================================================================
    est <- rowNames <- NULL

    # parameter estimates
    for (j in 1) {
      for (z in c(0,1)) {
        eval(parse(text = paste("est <- rbind(est, dat_david2$margF",
                                j, ".z", z, ".t0[1])", sep = "")))
        rowNames <- c(rowNames, paste(c(z, j), collapse = " "))
      }
    }
    row.names(est) <- rowNames
    var <- t(as.matrix(infCurves)) %*% as.matrix(infCurves)/n.data^2
    row.names(var) <- colnames(var) <- rowNames

    # output static interventions
    est <- 1 - est['1 1',]
    var <- var['1 1', '1 1']

    # if (all(dW == 1)) {
    #     est <- 1 - est['1 1',]
    #     var <- var['1 1', '1 1']
    # }else if (all(dW == 0)) {
    #     est <- 1 - est['0 1',]
    #     var <- var['0 1', '0 1']
    # }
    onestep_out_all[[tk_count]] <- list(est = est, var = var, meanIC = meanIC, ic = infCurves)
  }

  s_vec <- sapply(onestep_out_all, function(x) x$est)
  survival_df <- data.frame(s_vec, T.uniq)
  class(survival_df) <- 'surv_survtmle'

  return(list(survival_df = survival_df, onestep_out_all = onestep_out_all))
}



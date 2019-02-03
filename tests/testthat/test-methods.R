context("methods for the survival curve")
# set.seed(628957)
simulate_data <- function(n_sim = 2e2) {
  library(simcausal)
  D <- DAG.empty()
  D <- D +
    node("W1", distr = "rbinom", size = 1, prob = .5) +
    node("W", distr = "runif", min = 0, max = 1.5) +
    node("A", distr = "rbinom", size = 1, prob = .15 + .5 * as.numeric(W > .75)) +
    node("Trexp", distr = "rexp", rate = 1 + .7 * W^2 - .8 * A) +
    node("Cweib", distr = "rweibull", shape = 1 + .5 * W, scale = 75) +
    node("T", distr = "rconst", const = round(Trexp * 2)) +
    node("C", distr = "rconst", const = round(Cweib * 2)) +
    # Observed random variable (follow-up time):
    node("T.tilde", distr = "rconst", const = ifelse(T <= C, T, C)) +
    # Observed random variable (censoring indicator, 1 - failure event, 0 - censored):
    node("Delta", distr = "rconst", const = ifelse(T <= C, 1, 0))
  setD <- set.DAG(D)
  dat <- sim(setD, n = n_sim)
  # only grab ID, W's, A, T.tilde, Delta
  Wname <- grep("W", colnames(dat), value = TRUE)
  dat <- dat[, c("ID", Wname, "A", "T.tilde", "Delta")]
  # input: scalar q, W vector. computes for all W, the S(q|A,W)
  true_surv_one <- function(q, W, A = 1) sapply(W, function(w) {
      1 - pexp(q, rate = 1 + .7 * w^2 - .8 * A)
    })
  # input: vector q. mean(S(q|A,W)|A), average out W. loop over q
  true_surv <- function(q_grid, surv_fn, A) {
    W_grid <- seq(0, 1.5, .01)
    survout <- numeric()
    for (q in q_grid) survout <- c(survout, mean(surv_fn(q = q / 2, W = W_grid, A = A)))
    return(survout)
  }
  truth_surv <- function(q) true_surv(q_grid = q, surv_fn = true_surv_one, A = 1)
  truth_surv0 <- function(q) true_surv(q_grid = q, surv_fn = true_surv_one, A = 0)
  return(list(dat = dat, true_surv1 = truth_surv, true_surv0 = truth_surv0))
}
# simulation
n_sim <- 2e2
simulated <- simulate_data(n_sim = n_sim)
df <- simulated$dat
true_surv <- simulated$true_surv1

sl_lib_g <- c("SL.mean", "SL.glm")
sl_lib_censor <- c("SL.mean", "SL.glm")
sl_lib_failure <- c("SL.mean", "SL.glm", "SL.step.forward")
range(df$T.tilde)
df$T.tilde <- df$T.tilde + 1
k_grid <- 1:max(df$T.tilde)

sl_fit <- initial_sl_fit(
  ftime = df$T.tilde,
  ftype = df$Delta,
  trt = df$A,
  adjustVars = data.frame(df[, c("W", "W1")]),
  t_0 = max(df$T.tilde),
  SL.trt = sl_lib_g,
  SL.ctime = sl_lib_censor,
  SL.ftime = sl_lib_failure
)
sl_fit$density_failure_1$hazard_to_survival()
sl_fit$density_failure_0$hazard_to_survival()
sl_fit$density_failure_1$t <- k_grid
sl_fit$density_failure_0$t <- k_grid

# sl_density_failure_1_marginal <- sl_fit$density_failure_1$clone(deep = TRUE)
# sl_density_failure_0_marginal <- sl_fit$density_failure_0$clone(deep = TRUE)
# sl_density_failure_1_marginal$survival <- matrix(colMeans(sl_density_failure_1_marginal$survival), nrow = 1)
# sl_density_failure_0_marginal$survival <- matrix(colMeans(sl_density_failure_0_marginal$survival), nrow = 1)
test_that("sl_1 results should not be NA", {
  expect_true(all(!sapply(sl_fit$density_failure_1$survival, is.na)))
})
test_that("sl_0 results should not be NA", {
  expect_true(all(!sapply(sl_fit$density_failure_0$survival, is.na)))
})

################################################################################
# ipcw
ipcw_fit_1_all <- repeat_t_grid$new(
  method = ipcw,
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1
)$fit(k_grid = k_grid)
ipcw_fit_0_all <- repeat_t_grid$new(
  method = ipcw,
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_0,
  density_censor = sl_fit$density_censor_0,
  g1W = sl_fit$g1W,
  A_intervene = 0
)$fit(k_grid = k_grid)
ee_fit_1_all <- repeat_t_grid$new(
  method = ee,
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1
)$fit(k_grid = k_grid)
ee_fit_0_all <- repeat_t_grid$new(
  method = ee,
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_0,
  density_censor = sl_fit$density_censor_0,
  g1W = sl_fit$g1W,
  A_intervene = 0
)$fit(k_grid = k_grid)
ipcw_fit_1 <- survival_curve$new(t = k_grid, survival = ipcw_fit_1_all)
ipcw_fit_0 <- survival_curve$new(t = k_grid, survival = ipcw_fit_0_all)
ee_fit_1 <- survival_curve$new(t = k_grid, survival = ee_fit_1_all)
ee_fit_0 <- survival_curve$new(t = k_grid, survival = ee_fit_0_all)

test_that("ipcw_1 results should not be NA", {
  expect_true(all(!sapply(ipcw_fit_1$survival, is.na)))
})
test_that("ipcw_0 results should not be NA", {
  expect_true(all(!sapply(ipcw_fit_0$survival, is.na)))
})
test_that("ee_1 results should not be NA", {
  expect_true(all(!sapply(ee_fit_1$survival, is.na)))
})
test_that("ee_0 results should not be NA", {
  expect_true(all(!sapply(ee_fit_0$survival, is.na)))
})
################################################################################
# moss
moss_fit <- MOSS$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1,
  k_grid = k_grid
)
psi_moss_1 <- moss_fit$onestep_curve(
  epsilon = 1e-1 / n_sim,
  max_num_interation = 1e2,
  verbose = F
)
moss_fit <- MOSS$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_0,
  density_censor = sl_fit$density_censor_0,
  g1W = sl_fit$g1W,
  A_intervene = 0,
  k_grid = k_grid
)
psi_moss_0 <- moss_fit$onestep_curve(
  epsilon = 1e-1 / n_sim,
  max_num_interation = 1e2,
  verbose = F
)
moss_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_1)
moss_fit_0 <- survival_curve$new(t = k_grid, survival = psi_moss_0)
test_that("moss_1 results should not be NA", {
  expect_true(all(!sapply(moss_fit_1$survival, is.na)))
})
test_that("moss_0 results should not be NA", {
  expect_true(all(!sapply(moss_fit_0$survival, is.na)))
})

################################################################################
# moss hazard submodel
moss_hazard_fit <- MOSS_hazard$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1,
  k_grid = k_grid
)
psi_moss_hazard_1 <- moss_hazard_fit$iterate_onestep(epsilon = 1e-2, verbose = FALSE)
moss_hazard_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_1)
test_that("moss_hazard_1 results should not be NA", {
  expect_true(all(!sapply(moss_hazard_fit_1$survival, is.na)))
})

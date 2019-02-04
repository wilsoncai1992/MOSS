context("methods for the survival curve")
# set.seed(1234)
set.seed(11)
source("./simulate_data.R")
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
psi_moss_hazard_1 <- moss_hazard_fit$iterate_onestep(
  epsilon = 1e-2, verbose = FALSE
)
moss_hazard_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_1)
test_that("moss_hazard_1 results should not be NA", {
  expect_true(all(!sapply(moss_hazard_fit_1$survival, is.na)))
})

################################################################################
# moss difference curve

moss_hazard_ate_fit <- MOSS_hazard_ate$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  density_failure_0 = sl_fit$density_failure_0,
  density_censor_0 = sl_fit$density_censor_0,
  g1W = sl_fit$g1W,
  k_grid = k_grid
)
psi_moss_hazard_ate_1 <- moss_hazard_ate_fit$iterate_onestep(
  epsilon = 1e-2, max_num_interation = 5e1, verbose = FALSE
)
moss_hazard_ate_fit_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_ate_1)
test_that("moss_hazard_ate_1 results should not be NA", {
  expect_true(all(!sapply(moss_hazard_ate_fit_1$survival, is.na)))
})

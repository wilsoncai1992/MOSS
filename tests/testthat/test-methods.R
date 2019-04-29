context("methods for the survival curve")
# set.seed(11)
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
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  A = df$A,
  W = data.frame(df[, c("W", "W1")]),
  t_max = max(df$T.tilde),
  sl_treatment = sl_lib_g,
  sl_censoring = sl_lib_censor,
  sl_failure = sl_lib_failure
)
sl_fit$density_failure_1$hazard_to_survival()
sl_fit$density_failure_0$hazard_to_survival()
sl_fit$density_failure_1$t <- k_grid
sl_fit$density_failure_0$t <- k_grid

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
  max_num_interation = 1e1,
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
  max_num_interation = 1e1,
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
moss_hazard_l2 <- MOSS_hazard$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  A_intervene = 1,
  k_grid = k_grid
)
moss_hazard_l1 <- moss_hazard_l2$clone(deep = TRUE)
psi_moss_l2_1 <- moss_hazard_l2$iterate_onestep(
  method = "l2", epsilon = 1e-2, max_num_interation = 1e1, verbose = FALSE
)
moss_hazard_l2_1 <- survival_curve$new(t = k_grid, survival = psi_moss_l2_1)
test_that("MOSS l2 submodel results should not be NA", {
  expect_true(all(!sapply(moss_hazard_l2_1$survival, is.na)))
})
stoppingl2 <- mean(abs(moss_hazard_l2$compute_mean_eic(
  psi_n = psi_moss_l2_1,
  k_grid = moss_hazard_l2$q_best$t
)))
eic_fit <- eic$new(
  A = df$A,
  T_tilde = df$T.tilde,
  Delta = df$Delta,
  density_failure = sl_fit$density_failure_1,
  density_censor = sl_fit$density_censor_1,
  g1W = sl_fit$g1W,
  psi = colMeans(sl_fit$density_failure_1$survival),
  A_intervene = 1
)$all_t(k_grid = k_grid)
stopping0 <- mean(abs(colMeans(eic_fit)))
test_that("MOSS l2 submodel improves the stopping criteria", {
  expect_true(stoppingl2 < stopping0)
})

psi_moss_hazard_l1_1 <- moss_hazard_l1$iterate_onestep(
  method = "l1", epsilon = 1e-2, max_num_interation = 1e1, verbose = FALSE
)
moss_hazard_l1_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_l1_1)
test_that("MOSS l1 submodel results should not be NA", {
  expect_true(all(!sapply(moss_hazard_l1_1$survival, is.na)))
})

stoppingl1 <- mean(abs(moss_hazard_l1$compute_mean_eic(
  psi_n = psi_moss_hazard_l1_1,
  k_grid = moss_hazard_l1$q_best$t
)))
test_that("MOSS l1 submodel improves the stopping criteria", {
  expect_true(stoppingl1 < stopping0)
})

################################################################################
# moss difference curve

moss_hazard_ate_l2 <- MOSS_hazard_ate$new(
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
moss_hazard_ate_l1 <- moss_hazard_ate_l2$clone(deep = TRUE)
psi_moss_hazard_ate_l2_1 <- moss_hazard_ate_l2$iterate_onestep(
  method = "l2", epsilon = 1e-2, max_num_interation = 1e1, verbose = FALSE
)
moss_hazard_ate_fit_l2_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_ate_l2_1)
test_that("MOSS ATE l2 submodel should not be NA", {
  expect_true(all(!sapply(moss_hazard_ate_fit_l2_1$survival, is.na)))
})

psi_moss_hazard_ate_l1_1 <- moss_hazard_ate_l1$iterate_onestep(
  method = "l1", epsilon = 1e-2, max_num_interation = 1e1, verbose = FALSE
)
moss_hazard_ate_fit_l1_1 <- survival_curve$new(t = k_grid, survival = psi_moss_hazard_ate_l1_1)
test_that("MOSS ATE l1 submodel should not be NA", {
  expect_true(all(!sapply(moss_hazard_ate_fit_l1_1$survival, is.na)))
})

context("methods for the survival curve")
# set.seed(628957)
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

sl_density_failure_1_marginal <- sl_fit$density_failure_1$clone(deep = TRUE)
sl_density_failure_0_marginal <- sl_fit$density_failure_0$clone(deep = TRUE)
sl_density_failure_1_marginal$survival <- matrix(
  colMeans(sl_density_failure_1_marginal$survival), nrow = 1
)
sl_density_failure_0_marginal$survival <- matrix(
  colMeans(sl_density_failure_0_marginal$survival), nrow = 1
)

sl_fit$density_failure_1$survival_to_hazard()
sl_fit$density_failure_1$hazard_to_pdf()
sl_fit$density_failure_1$pdf_to_survival()
sl_fit$density_failure_1$pdf_to_hazard()

context("constrained enet regression is working")
library(Rsolnp)
library(glmnet)
n <- 1e4
X <- matrix(rnorm(13 * n), ncol = 13)
beta <- c(c(10, 5, 1), rep(0, 10))
pY <- MOSS:::expit(X %*% beta)
Y <- rbinom(n = n, size = 1, prob = pY)

glmnet_fit <- glmnet::glmnet(
  x = X,
  y = Y,
  family = "binomial",
  alpha = 0,
  lambda = 1e-5,
  intercept = FALSE,
  standardize = FALSE
)
sol1 <- as.vector(glmnet_fit$beta)

glmnet_cv_fit <- glmnet::cv.glmnet(
  x = X,
  y = Y,
  family = "binomial",
  alpha = 0,
  intercept = FALSE,
  standardize = FALSE,
  lambda = 10 ^ seq(-9, -3, length.out = 1e2)
)
sol1_2 <- coef(glmnet_cv_fit)[-1]

l2_norm_max <- 16
sol2 <- fit_enet_constrained(
  Y = Y, X = X, beta_init = rep(1, 13), norm_max = l2_norm_max, type = "l2"
)


test_that("classic glmnet is working", {
  expect_true(sum(abs(beta - sol1)) < 2)
})
test_that("constrained ridge regression is working", {
  expect_true(sum(abs(beta - sol2)) < 2)
})
test_that("constrained ridge regression is close to glmnet", {
  expect_true(sum(abs(sol1 - sol2)) < 1)
})
test_that("constrained ridge regression is obeying constraint", {
  expect_true(sqrt(sum(sol2 ^ 2)) <= l2_norm_max + 1e-3)
})

l1_norm_max <- 16
sol_lasso2 <- fit_enet_constrained(
  Y = Y, X = X, beta_init = rep(1, 13), norm_max = l1_norm_max, type = "l1"
)
test_that("constrained l1 regression is close to glmnet", {
  expect_true(sum(abs(beta - sol_lasso2)) < 2)
})

test_that("constrained l1 regression is obeying constraint", {
  expect_true(sum(abs(sol_lasso2)) <= l1_norm_max + 1e-3)
})

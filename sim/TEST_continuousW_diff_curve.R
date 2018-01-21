# simulate data
library(simcausal)
D <- DAG.empty()
D <- D +
  node("W", distr = "runif", min = 0, max = 1) +
  node("A", distr = "rbinom", size = 1, prob = .2 + .5*W) +
  # node("A", distr = "rbinom", size = 1, prob = .5) +
  node("Trexp", distr = "rexp", rate = 1 + .5*W - .5*A) +
  node("Cweib", distr = "rweibull", shape = .7 - .2*W, scale = 1) +
  node("T", distr = "rconst", const = round(Trexp*10,0)) +
  node("C", distr = "rconst", const = round(Cweib*10, 0)) +
  # node("C", distr = "rconst", const = 999) +
  node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
  node("Delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
setD <- set.DAG(D)

# Simulate the data from the above data generating distribution:
# dat <- sim(setD, n=1e2, rndseed = 12345)
dat <- sim(setD, n=1e3, rndseed = 12345)
# dat <- sim(setD, n=1e4, rndseed = 12345)
head(dat)

# subset into observed dataset
library(dplyr)
# only grab ID, W's, A, T.tilde, Delta
Wname <- grep('W', colnames(dat), value = TRUE)
dat <- dat[,c('ID', Wname, 'A', "T.tilde", "Delta")]
head(dat)
# ---------------------------------------------------------------------------------------
# KM
library(survival)
n.data <- nrow(dat)
km.fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = dat)
plot(km.fit, col=c('blue','red'), lty = 2,
main = paste('n=', n.data, '\n # of nuisance covariate = 1'),
xlab = 'Time')

library(dplyr)
true_func <- function(x,A,W){
  1 - pexp(x/10, rate = 1+.5*W-.5*A)
}
plot_one_arm <- function(A, ...) {
  meshgrid <- expand.grid(seq(0,50,.1), A, seq(0,1,.1))
  p_hat <- apply(meshgrid, 1, function(x) true_func(x[1],x[2],x[3]))
  df_temp <- data.frame(meshgrid, p_hat)
  colnames(df_temp)[1] <- 'x'
  df_temp <- df_temp %>% group_by(x) %>% summarise(p = mean(p_hat)) %>% as.data.frame
  lines(df_temp[,1], df_temp[,2], type="l", cex=0.2, ...)
}
plot_one_arm(A = 1, col = 'red')
plot_one_arm(A = 0, col = 'blue')
# ---------------------------------------------------------------------------------------
library(MOSS)
# R6
onestepfit = MOSS_difference$new(dat,
  verbose = TRUE, epsilon.step = 1e-2, max.iter = 1e2)
  # verbose = TRUE, epsilon.step = 1e-2, max.iter = 2e2)
onestepfit$initial_fit()
onestepfit$onestep_diff_curve()
onestepfit$print()
onestepfit$plot_CI_pointwise()
onestepfit$plot_CI_pointwise(ylim = c(-.5,.5))

onestepfit$compute_CI_simultaneous()
onestepfit$plot_CI_simultaneous(ylim = c(-.5,.5))
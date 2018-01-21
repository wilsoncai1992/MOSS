# ====================================================================================================
# simulate data
# ====================================================================================================
library(simcausal)
D <- DAG.empty()

D <- D +
  node("W", distr = "rbinom", size = 1, prob = .5) +
  node("A", distr = "rbinom", size = 1, prob = .3 + .3*W) +
  # node("A", distr = "rbinom", size = 1, prob = .5) +
  node("Trexp", distr = "rexp", rate = 1 + .5*W - .5*A) +
  # node("Trexp", distr = "rexp", rate = 1 + 2*W - .5*A) +
  node("Cweib", distr = "rweibull", shape = .7 - .2*W, scale = 1) +
  node("T", distr = "rconst", const = round(Trexp*10,0)) +
  # node("C", distr = "rconst", const = round(Cweib*10, 0)) +
  node("C", distr = "rconst", const = 999) +
  node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
  node("Delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
setD <- set.DAG(D)

# Simulate the data from the above data generating distribution:
dat <- sim(setD, n=1e2, rndseed = 12345)
# dat <- sim(setD, n=1e3, rndseed = 12345)
# dat <- sim(setD, n=1e4, rndseed = 12345)
head(dat)

# subset into observed dataset
library(dplyr)
# only grab ID, W's, A, T.tilde, Delta
Wname <- grep('W', colnames(dat), value = TRUE)
dat <- dat[,c('ID', Wname, 'A', "T.tilde", "Delta")]
head(dat)

# plot KM
library(survival)
n.data <- nrow(dat)
km.fit <- survfit(Surv(time = T.tilde, event = Delta) ~ A, data = dat)
plot(km.fit, col=c('blue','red'), lty = 2,
main = paste('n=', n.data, '\n # of nuisance covariate = 1'),
xlab = 'Time')

q <- seq(0,3,.01)
# truesurvExp1 <- 1 - pexp(q, rate = 1.5)
# truesurvExp2 <- 1 - pexp(q, rate = 1)
truesurvExp1 <- 1 - pexp(q, rate = 1)
truesurvExp2 <- 1 - pexp(q, rate = .5)
truesurvExp <- (truesurvExp1 + truesurvExp2)/2
lines(round(q*10,0), truesurvExp, type="l", cex=0.2, col = 'red')

truesurvExp1 <- 1 - pexp(q, rate = 1.5)
truesurvExp2 <- 1 - pexp(q, rate = 1)
# truesurvExp1 <- 1 - pexp(q, rate = 2)
# truesurvExp2 <- 1 - pexp(q, rate = 1.5)
# truesurvExp1 <- 1 - pexp(q, rate = 3)
# truesurvExp2 <- 1 - pexp(q, rate = 2.5)
truesurvExp <- (truesurvExp1 + truesurvExp2)/2
lines(round(q*10,0), truesurvExp, type="l", cex=0.2, col = 'blue')



# R6
# onestepfit = MOSS$new(dat, dW = 1, verbose = TRUE, epsilon.step = 1e-3, max.iter = 1e2)
onestepfit = MOSS$new(dat, dW = 1,
  # verbose = TRUE, epsilon.step = 1e-3, max.iter = 5e2)
  verbose = TRUE, epsilon.step = 1e-1, max.iter = 5e2)
  # verbose = TRUE, epsilon.step = 5e-4, max.iter = 5e2)
  # verbose = TRUE, epsilon.step = 1e-4, max.iter = 5e2)
# onestepfit$initial_fit(g.SL.Lib = c("SL.glm", "SL.step", "SL.glm.interaction"),
onestepfit$initial_fit(g.SL.Lib = c("SL.mean","SL.glm", 'SL.gam'),
                       Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam"),
                       ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam"))
onestepfit$transform_failure_hazard_to_survival()
onestepfit$transform_failure_hazard_to_pdf()
onestepfit$compute_EIC()


iter_count <- 0
stopping_prev <- Inf
stopping_history <- numeric()
all_loglikeli <- numeric()

stopping <- onestepfit$compute_stopping()
while ((stopping >= onestepfit$tol) & (iter_count <= onestepfit$max.iter)) {
# while ((stopping >= onestepfit$tol) & (iter_count <= onestepfit$max.iter) & ((stopping_prev - stopping) >= max(-onestepfit$tol, -1e-5))) {
  print(stopping)
  # onestepfit$onestep_curve_update()
  # onestepfit$onestep_curve_update_pooled()
  onestepfit$onestep_curve_update_mat()
  onestepfit$compute_EIC()
  iter_count <- iter_count + 1
  stopping_history[iter_count] <- stopping
  stopping_prev <- stopping
  stopping <- onestepfit$compute_stopping()

  onestepfit$compute_Psi()
  if (iter_count %% 10 == 0) onestepfit$plot_onestep_curve(add = TRUE)
  # if (iter_count %% 10 == 0) plot(onestepfit$Pn.D1.t); abline(h = 0)
}

if (iter_count == onestepfit$max.iter) {
  warning('Max Iter count reached, stop iteration.')
}

onestepfit$compute_Psi()
onestepfit$Psi.hat

onestepfit$plot_CI_pointwise(add = TRUE)

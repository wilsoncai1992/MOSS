# ====================================================================================================
# simulate data
# ====================================================================================================
library(simcausal)
D <- DAG.empty()

D <- D +
  node("W1", distr = "rbinom", size = 1, prob = .5) +
  node("W", distr = "rbinom", size = 1, prob = .5) +
  node("A", distr = "rbinom", size = 1, prob = .15 + .5*W) +
  node("Trexp", distr = "rexp", rate = 1 + .5*W - .5*A) +
  node("Cweib", distr = "rweibull", shape = .7 - .2*W, scale = 1) +
  node("T", distr = "rconst", const = round(Trexp*100,0)) +
  node("C", distr = "rconst", const = round(Cweib*100, 0)) +
  node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
  node("Delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
setD <- set.DAG(D)

# Simulate the data from the above data generating distribution:
dat <- sim(setD, n=3e2)
head(dat)

# subset into observed dataset
library(dplyr)
# only grab ID, W's, A, T.tilde, Delta
Wname <- grep('W', colnames(dat), value = TRUE)
dat <- dat[,c('ID', Wname, 'A', "T.tilde", "Delta")]
head(dat)

# R6
onestepfit = MOSS$new(dat, dW = 1, verbose = TRUE, epsilon.step = 1e-3, max.iter = 1e3)
onestepfit$display()

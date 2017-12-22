# simulate data
set.seed(341796)
n <- 1e3
t_0 <- 100
W <- data.frame(W1 = runif(n), W2 = rbinom(n, 1, 0.5))
A <- rbinom(n, 1, 0.5)
T <- rgeom(n,plogis(-4 + W$W1 * W$W2 - A)) + 1
C <- rgeom(n, plogis(-6 + W$W1)) + 1
ftime <- pmin(T, C)
ftype <- as.numeric(ftime == T)

# load the package
library(survtmle)
#> survtmle: Targeted Learning for Survival Analysis
#> Version: 1.0.2.1

# apply survtmle for estimation
fit <- survtmle(ftime = ftime, ftype = ftype,
                trt = A, adjustVars = W,
                SL.trt = c('SL.glm'),
                # glm.trt = "1",
                SL.ftime = c('SL.glm'),
                # glm.ftime = "I(W1*W2) + trt + t",
                SL.ctime = c('SL.glm'),
                # glm.ctime = "W1 + t",
                method = "hazard",
                t0 = t_0)

# extract cumulative incidence at each timepoint
tpfit <- timepoints(fit, times = seq_len(t_0))

# examine output object produced by the timepoints function
tpfit
# examine plot of cumulative incidences
plot(tpfit)

# extract two conterfactual survival curve
len_groups <- as.numeric(unique(lapply(lapply(tpfit, FUN = `[[`,
    "est"), FUN = length)))
names_groups <- unique(lapply(lapply(tpfit, FUN = `[[`, "est"),
    FUN = rownames))[[1]]
est_only <- t(matrix(unlist(lapply(tpfit, FUN = `[[`, "est")),
    ncol = len_groups, byrow = TRUE))
est_only <- as.data.frame(est_only)
rownames(est_only) <- names_groups
colnames(est_only) <- paste0("t", seq_len(ncol(est_only)))

s_0 <- 1 - as.numeric(est_only[1,])
s_1 <- 1 - as.numeric(est_only[2,])

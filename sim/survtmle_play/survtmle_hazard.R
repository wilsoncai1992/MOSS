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
#
# fit hazard
# ---------------------------------------------------------------------------------------
adjustVars <- data.frame(W)
trt <- A
ftypeOfInterest <- unique(ftype)
trtOfInterest <- 0:1
SL.trt <- c('SL.glm')
gtol <- 1e-3
SL.ctime <- c('SL.glm')
SL.ftime <- c('SL.glm')

n <- length(ftime)
id <- seq_len(n)
dat <- data.frame(id = id, ftime = ftime, ftype = ftype, trt = trt)
if(!is.null(adjustVars)) dat <- cbind(dat, adjustVars)

nJ <- length(ftypeOfInterest)
allJ <- sort(unique(ftype[ftype != 0]))
ofInterestJ <- sort(ftypeOfInterest)

# calculate number of groups
ntrt <- length(trtOfInterest)
uniqtrt <- sort(trtOfInterest)

# estimate trt probabilities
trtOut <- survtmle:::estimateTreatment(dat = dat,
                            ntrt = ntrt,
                            uniqtrt = uniqtrt,
                            adjustVars = adjustVars,
                            SL.trt = SL.trt,
                            # glm.trt = glm.trt,
                            returnModels = TRUE,
                            gtol = gtol)
dat <- trtOut$dat
trtMod <- trtOut$trtMod

# make long version of data sets needed for estimation and prediction
dataList <- survtmle:::makeDataList(dat = dat, J = allJ, ntrt = ntrt, uniqtrt = uniqtrt,
                         t0 = t_0, bounds = NULL)

# estimate censoring
censOut <- survtmle:::estimateCensoring(dataList = dataList,
                             ntrt = ntrt,
                             uniqtrt = uniqtrt,
                             t0 = t_0,
                             verbose = FALSE,
                             adjustVars = adjustVars,
                             SL.ctime = SL.ctime,
                             # glm.ctime = glm.ctime,
                             glm.family = 'binomial',
                             returnModels = TRUE,
                             gtol = gtol)
dataList <- censOut$dataList
ctimeMod <- censOut$ctimeMod

# estimate cause specific hazards
estOut <- survtmle:::estimateHazards(dataList = dataList,
                          J = allJ,
                          verbose = FALSE,
                          bounds = NULL,
                          adjustVars = adjustVars,
                          SL.ftime = SL.ftime,
                          # glm.ftime = glm.ftime,
                          glm.family = 'binomial',
                          returnModels = TRUE)
dataList <- estOut$dataList
ftimeMod <- estOut$ftimeMod
# check for convergence
suppressWarnings(
  if(all(dataList[[1]] == "convergence failure")) {
    return("estimation convergence failure")
  }
)

# extract g
g_1 <- dat$g_1
g_0 <- dat$g_0

# extract hazard
d1 <- dataList$`1`
d0 <- dataList$`0`

haz1 <- d1[,c('id', 't', 'Q1Haz')]
haz1 <- tidyr::spread(haz1, t, Q1Haz)
haz1 <- haz1[,-1] # remove the id column

haz0 <- d0[,c('id', 't', 'Q1Haz')]
haz0 <- tidyr::spread(haz0, t, Q1Haz)
haz0 <- haz0[,-1] # remove the id column

# extract S_{Ac}
S_Ac_1 <- d1[,c('id', 't', 'G_dC')]
S_Ac_1 <- tidyr::spread(S_Ac_1, t, G_dC)
S_Ac_1 <- S_Ac_1[,-1] # remove the id column

S_Ac_0 <- d0[,c('id', 't', 'G_dC')]
S_Ac_0 <- tidyr::spread(S_Ac_0, t, G_dC)
S_Ac_0 <- S_Ac_0[,-1] # remove the id column


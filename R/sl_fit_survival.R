#' super learner fit for failure and censoring event
#'
#' using survtmle package
#'
#' @param ftime vector of last follow up time
#' @param ftype vector of censoring indicator
#' @param trt vector of treatment
#' @param adjustVars data.frame of baseline covariates
#' @param t_0 the maximum time to estimate the survival probabilities
#' @param trtOfInterest the intervention of interest
#' @param SL.trt SuperLearner library for propensity score
#' @param gtol treshold for the fitted propensity scores
#' @param SL.ctime SuperLearner library for censoring event hazard
#' @param SL.ftime SuperLearner library for failure event hazard
#'
#' @importFrom SuperLearner SuperLearner
#' @importFrom survtmle estimateTreatment makeDataList estimateCensoring estimateHazards
#' @export
initial_sl_fit <- function(
                           ftime,
                           ftype,
                           trt,
                           adjustVars,
                           t_0,
                           trtOfInterest = 0:1,
                           SL.trt = c("SL.glm"),
                           gtol = 1e-3,
                           SL.ctime = c("SL.glm"),
                           SL.ftime = c("SL.glm")) {
  adjustVars <- data.frame(adjustVars)
  ftypeOfInterest <- unique(ftype)
  n <- length(ftime)
  id <- seq_len(n)
  dat <- data.frame(id = id, ftime = ftime, ftype = ftype, trt = trt)
  if (!is.null(adjustVars)) dat <- cbind(dat, adjustVars)

  nJ <- length(ftypeOfInterest)
  allJ <- sort(unique(ftype[ftype != 0]))
  ofInterestJ <- sort(ftypeOfInterest)

  # calculate number of groups
  ntrt <- length(trtOfInterest)
  uniqtrt <- sort(trtOfInterest)

  # estimate trt probabilities
  trtOut <- survtmle::estimateTreatment(
    dat = dat,
    ntrt = ntrt,
    uniqtrt = uniqtrt,
    adjustVars = adjustVars,
    SL.trt = SL.trt,
    returnModels = TRUE,
    gtol = gtol
  )
  dat <- trtOut$dat
  trtMod <- trtOut$trtMod

  # make long version of data sets needed for estimation and prediction
  dataList <- survtmle::makeDataList(
    dat = dat, J = allJ, ntrt = ntrt, uniqtrt = uniqtrt, t0 = t_0, bounds = NULL
  )
  # estimate censoring
  # when there is almost no censoring, the classification will fail;
  # we manually input the conditional survival for the censoring
  censOut <- tryCatch({
    survtmle::estimateCensoring(
      dataList = dataList,
      ntrt = ntrt,
      uniqtrt = uniqtrt,
      t0 = t_0,
      verbose = FALSE,
      adjustVars = adjustVars,
      SL.ctime = SL.ctime,
      glm.family = "binomial",
      returnModels = TRUE,
      gtol = gtol
    )
  },
  error = function(cond) {
    message("censoring sl error")
    NULL
  }
  )
  if (is.null(censOut)) {
    censOut <- list()
    censOut$dataList <- dataList
    censOut$dataList$obs[, "G_dC"] <- 1
    censOut$dataList$'0'[, "G_dC"] <- 1
    censOut$dataList$'1'[, "G_dC"] <- 1
    is_sl_censoring_converge <- FALSE
    dataList <- censOut$dataList
  } else {
    dataList <- censOut$dataList
    ctimeMod <- censOut$ctimeMod
    is_sl_censoring_converge <- TRUE
  }

  # estimate cause specific hazards
  estOut <- survtmle::estimateHazards(
    dataList = dataList,
    J = allJ,
    verbose = FALSE,
    bounds = NULL,
    adjustVars = adjustVars,
    SL.ftime = SL.ftime,
    glm.family = "binomial",
    returnModels = TRUE
  )
  dataList <- estOut$dataList
  ftimeMod <- estOut$ftimeMod
  # check for convergence
  suppressWarnings(
    if (all(dataList[[1]] == "convergence failure")) {
      return("estimation convergence failure")
    }
  )

  # extract g
  g_1 <- dat$g_1
  g_0 <- dat$g_0

  # extract hazard
  d1 <- dataList$`1`
  d0 <- dataList$`0`

  haz1 <- d1[, c("id", "t", "Q1Haz")]
  haz1 <- tidyr::spread(haz1, t, Q1Haz)
  haz1$id <- NULL # remove the id column

  haz0 <- d0[, c("id", "t", "Q1Haz")]
  haz0 <- tidyr::spread(haz0, t, Q1Haz)
  haz0$id <- NULL # remove the id column

  # extract S_{Ac}
  S_Ac_1 <- d1[, c("id", "t", "G_dC")]
  S_Ac_1 <- tidyr::spread(S_Ac_1, t, G_dC)
  S_Ac_1 <- S_Ac_1[, -1] # remove the id column

  S_Ac_0 <- d0[, c("id", "t", "G_dC")]
  S_Ac_0 <- tidyr::spread(S_Ac_0, t, G_dC)
  S_Ac_0 <- S_Ac_0[, -1] # remove the id column

  density_failure_1 <- survival_curve$new(
    t = seq(range(ftime)[1], range(ftime)[2]), hazard = haz1
  )
  density_failure_0 <- survival_curve$new(
    t = seq(range(ftime)[1], range(ftime)[2]), hazard = haz0
  )
  density_censor_1 <- survival_curve$new(
    t = seq(range(ftime)[1], range(ftime)[2]), survival = S_Ac_1
  )
  density_censor_0 <- survival_curve$new(
    t = seq(range(ftime)[1], range(ftime)[2]), survival = S_Ac_0
  )
  return(list(
    density_failure_1 = density_failure_1,
    density_failure_0 = density_failure_0,
    density_censor_1 = density_censor_1,
    density_censor_0 = density_censor_0,
    g1W = g_1[, 1]
  ))
}

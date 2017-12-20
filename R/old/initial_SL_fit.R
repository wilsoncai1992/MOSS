#' wrapper for conditonal hazard regression SuperLearner
#'
#' use T, A, C, W data format as input
#'
#' @param dat data.frame, with col: T, A, C, W
#' @param T.uniq vector of unique event times
#' @param ht.SL.Lib library for SuperLearner
#'
#' @return h.hat.t matrix, nrow = n.sample, ncol = length(T.uniq);
#'          each row is the predicted conditional hazard h(T=t | T>=t, A, W) for that subject
#' @export
#'
#' @examples
#' # TO DO
#' @import dplyr
#' @importFrom tidyr spread
estimate_hazard_SL <- function(dat,
                               T.uniq,
                               ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth")) {
  # transform original data into SL-friendly format
  dat_david <- dat
  dat_david <- rename(dat_david, Z = A)

  if ('ID' %in% toupper(colnames(dat_david))) {
    # if there are already id in the dataset
    dat_david <- rename(dat_david, id = ID)
  }else{
    # if no id exist
    # create 'id' on our own
    dat_david$id <- 1:nrow(dat_david)
  }

  # censoring
  dat_david <- rename(dat_david, Delta.J = Delta)

  # remove all other useless columns
  baseline_name <- grep('W', colnames(dat_david), value = TRUE)
  keeps <- c("id", baseline_name, 'T.tilde', 'Delta.J', 'Z')
  dat_david <- dat_david[,keeps]

  # first fit SL on the max time point
  T.it <- max(T.uniq)
  # =======================================================================================
  # dataList
  # =======================================================================================
  datalist <- makeDataList(dat = dat_david,
                           J = 1, # one kind of failure
                           nZ = 2, # one kind of treatment
                           Z = c(0,1),
                           t0 = T.it, # time to predict on
                           bounds=NULL)
  # ----------------------------------------------------------------------------------------------------
  # perform hazard SL for the maximum time point
  Haz_hat <- estimateHazards(dataList = datalist, J=1, verbose=FALSE, strata=NULL,
                             adjustVars = dat_david[,baseline_name],
                             SLlibrary.event = ht.SL.Lib,
                             glmFormula.event = NULL,
                             bounds = NULL)
  # Haz_hat[[1]][Haz_hat[[1]]$t == T.it,'Q1Haz']
  # mean(Haz_hat[[1]][Haz_hat[[1]]$t == T.it,'Q1Haz'])

  # turn into wide format

  out_haz <- Haz_hat[[3]]
  out_haz <- out_haz[,c('id', 't', 'Q1Haz')]
  out_haz <- tidyr::spread(out_haz, t, Q1Haz)
  # the colname number correspond to h_{T>t-1 | T>=t-1}
  rownames(out_haz) <- out_haz$id

  # turn NA entries (after failure) into zero hazard
  out_haz[is.na(out_haz)] <- 0

  # remove the id column
  out_haz_2 <- out_haz[,-1]

  # subset the columns for those only in T.uniq
  out_haz <- out_haz_2[,T.uniq]

  return(list(out_haz = out_haz,
              out_haz_full = out_haz_2))
}



#' wrapper for censoring regression SuperLearner
#'
#' @param dat data.frame, with col: T, A, W, Delta
#' @param T.uniq the unique time points of failure
#' @param Delta.SL.Lib library for SuperLearner
#'
#' @return h.hat.t matrix, nrow = n.sample, ncol = length(T.uniq);
#' each row is the predicted conditional hazard h(T=t | T>=t, A, W) for that subject
#' @export
#'
#' @examples
#' # TO DO
#' @import SuperLearner
#' @import survtmle2
#' @import dplyr
#' @importFrom tidyr spread
estimate_censoring_SL <- function(dat,
                                  T.uniq,
                                  Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth")) {
  # transform original data into SL-friendly format
  dat_david <- dat

  dat_david <- rename(dat_david, ftime = T.tilde)
  dat_david <- rename(dat_david, trt = A)

  if ('ID' %in% toupper(colnames(dat_david))) {
    # if there are already id in the dataset
    dat_david <- rename(dat_david, id = ID)
  }else{
    # if no id exist
    # create 'id' on our own
    dat_david$id <- 1:nrow(dat_david)
  }

  # censoring
  if ('Delta' %in% colnames(dat_david)) dat_david <- rename(dat_david, ftype = Delta)

  # remove all other useless columns
  baseline_name <- grep('W', colnames(dat_david), value = TRUE)
  keeps <- c("id", baseline_name, 'ftime', 'ftype', 'trt')
  dat_david <- dat_david[, keeps]

  T.uniq <- unique(sort(dat_david$ftime))
  T.max <- max(T.uniq)

  adjustVars <- dat_david[,baseline_name]
  # ====================================================================================================
  # dataList
  # ====================================================================================================
  datalist <- survtmle2:::makeDataList(dat = dat_david,
                                      J = 1, # one kind of failure
                                      ntrt = 2, # one kind of treatment
                                      uniqtrt = c(0,1),
                                      t0 = T.max, # time to predict on
                                      bounds=NULL)
  # ====================================================================================================
  # estimate g_2 (censoring)
  # perform censoring SL for the maximum time point
  # ====================================================================================================
  g2_hat <- survtmle2:::estimateCensoring(dataList = datalist, adjustVars = adjustVars,
                                         t0 = T.max,
                                         ntrt = 2, # one kind of treatment
                                         uniqtrt = c(0,1),
                                         SL.ctime = Delta.SL.Lib,
                                         returnModels = FALSE,verbose = FALSE)
  dataList2 <- g2_hat$dataList
  # ----------------------------------------------------------------------------------------------------
  # turn into wide format
  out_censor <- dataList2$`1`
  out_censor <- out_censor[,c('id', 't', 'G_dC')]
  out_censor <- tidyr::spread(out_censor, t, G_dC)
  # the colname number correspond to h_{T>t-1 | T>=t-1}
  rownames(out_censor) <- out_censor$id

  # turn NA entries (after failure) into zero hazard
  out_censor[is.na(out_censor)] <- 0

  # remove the id column
  out_censor_2 <- out_censor[,-1]

  # subset the columns for those only in T.uniq
  out_censor <- out_censor_2[,T.uniq]

  return(list(out_censor = out_censor,
              out_censor_full = out_censor_2))
}


#' generate data.frame to be used in survtmle package
#'
#' @param dat a data.frame with columns named:
#'         id = unique subject id
#'         T.tilde = (possibly censored) failure time
#'         Delta.J= (possible censored) failure type -- can set all to 1 if no censoring
#'         Z = binary treatment indicator
#'         W = other columns for baseline variables
#' @param J the unique positive values of dat$Delta.J (set J=1 if only one type of failure)
#' @param nZ set to length of unique(dat$Z)
#' @param Z the treatment values that are of interest (e.g., Z=c(0,1))
#' @param t0 the time point you eventually want predictions at
#' @param bounds just leave NULL
#'
#' @return
#' @export
#'
#' @examples
#' # TO DO
makeDataList <- function(dat, J, nZ, Z, t0, bounds=NULL){
  n <- nrow(dat)
  dataList <- vector(mode="list",length=nZ+1)

  # first element used for estimation
  dataList[[1]] <- dat[rep(1:nrow(dat),dat$T.tilde),]
  for(j in J){
    eval(parse(text=paste("dataList[[1]]$N",j," <- 0",sep="")))
    eval(parse(text=paste("dataList[[1]]$N",j,"[cumsum(dat$T.tilde)] <- as.numeric(dat$Delta.J==j)",sep="")))
  }
  dataList[[1]]$C <- 0
  dataList[[1]]$C[cumsum(dat$T.tilde)] <- as.numeric(dat$Delta.J==0)

  n.row.ii <- nrow(dataList[[1]])
  row.names(dataList[[1]])[row.names(dataList[[1]]) %in% paste(row.names(dat))] <- paste(row.names(dat),".0",sep="")
  dataList[[1]]$t <- as.numeric(paste(unlist(strsplit(row.names(dataList[[1]]),".",fixed=T))[seq(2,n.row.ii*2,2)]))+1

  if(!is.null(bounds)){
    boundFormat <- data.frame(t=bounds$t)
    for(j in J){
      if(paste("l",j,sep="") %in% names(bounds)){
        eval(parse(text=paste("boundFormat$l",j," <- bounds$l",j,sep="")))
      }else{
        eval(parse(text=paste("boundFormat$l",j," <- 0",sep="")))
      }
      if(paste("u",j,sep="") %in% names(bounds)){
        eval(parse(text=paste("boundFormat$u",j," <- bounds$u",j,sep="")))
      }else{
        eval(parse(text=paste("boundFormat$u",j," <- 1",sep="")))
      }
    }
    suppressMessages(
      dataList[[1]] <- join(x=dataList[[1]],y=boundFormat,type="left")
    )
  }else{
    for(j in J){
      eval(parse(text=paste("dataList[[1]]$l",j," <- 0",sep="")))
      eval(parse(text=paste("dataList[[1]]$u",j," <- 1",sep="")))
    }
  }

  # subsequent elements used for prediction
  for(i in 1:nZ){
    dataList[[i+1]] <- dat[sort(rep(1:nrow(dat),t0)),]
    dataList[[i+1]]$t <- rep(1:t0,n)
    for(j in J){
      typejEvents <- dat$id[which(dat$Delta.J==j)]
      eval(parse(text=paste("dataList[[i+1]]$N",j," <- 0",sep="")))
      eval(parse(text=paste("dataList[[i+1]]$N",j,"[dataList[[i+1]]$id %in% typejEvents &  dataList[[i+1]]$t >= dataList[[i+1]]$T.tilde] <- 1",sep="")))
    }
    censEvents <- dat$id[which(dat$Delta.J==0)]
    dataList[[i+1]]$C <- 0
    dataList[[i+1]]$C[dataList[[i+1]]$id %in% censEvents & dataList[[i+1]]$t >= dataList[[i+1]]$T.tilde] <- 1
    dataList[[i+1]]$Z <- Z[i]
    dataList[[i+1]]$T.tilde <- t0

    if(!is.null(bounds)){
      suppressMessages(
        dataList[[i+1]] <- join(x=dataList[[i+1]],y=boundFormat,type="left")
      )
    }else{
      for(j in J){
        eval(parse(text=paste("dataList[[",i,"+1]]$l",j," <- 0",sep="")))
        eval(parse(text=paste("dataList[[",i,"+1]]$u",j," <- 1",sep="")))
      }
    }
  }
  names(dataList) <- c("obs",Z)
  return(dataList)
}

#' Conditional Hazard estimation in survtmle package
#' From survtmle (David Benkeser, 2016)
#'
#' @param dataList output of makeDataList
#' @param J same as above
#' @param verbose print messages as SL runs?
#' @param strata would leave as NULL, dont remember exactly how it is coded
#' @param adjustVars data.frame of predictors for SL, dont know why I have it coded this way...
#'                 just put in e.g. dat[,c("W1","W2")]. by default it will add the Z column
#'                 and t column in dataList[[1]] when it fits the SL
#' @param SLlibrary.event Super Learner library
#' @param glmFormula.event a glm formula if you don't want to use SL
#' @param superLearnerSummary doesn't do anything, just leave blank
#' @param bounds NULL
#'
#' @return  added column for the initial data.frame
#' QjHaz <- result, for each subject, h_t from time 0 to T-1; organized in T_i rows
#' QjPseudoHaz <-
#' returns list of length = length(dataList) with added columns name QjPseudoHaz and
#' QjHaz for all j in J. Ignore the PseudoHaz variables -- they're used when there are
#' multiple event types; if there's only one type of event they'll == Haz values anyway

#' @export
#'
#' @examples
#' # TO DO
estimateHazards <- function(dataList,
                            J,
                            verbose,
                            strata=NULL,
                            adjustVars,
                            SLlibrary.event,
                            glmFormula.event,
                            superLearnerSummary,
                            bounds){
  # check for missing inputs
  if(is.null(SLlibrary.event) & is.null(glmFormula.event) & is.null(strata)){
    warning("Super Learner library, glm formula, and strata for events not specified. Proceeding
            with empirical estimates")
    glmFormula.event <- "Z*factor(t)"
  }

  if(is.null(SLlibrary.event) & is.null(strata)){
    if(is.null(bounds)){
      for(j in J){
        # formula
        Qj.form <- sprintf("%s ~ %s", paste("N",j,sep=""), glmFormula.event)

        # add up all events less than current j to see who to include in regression
        NlessthanJ <- rep(0, nrow(dataList[[1]]))
        for(i in J[J<j]){
          eval(parse(text=paste("NlessthanJ <- NlessthanJ + dataList[[1]]$N",i,sep="")))
        }

        # fit glm
        Qj.mod <- glm(as.formula(Qj.form), data=dataList[[1]][NlessthanJ==0,], family="binomial")

        # get predictions back
        dataList <- lapply(dataList, function(x,j){
          eval(parse(text=paste("x$Q",j,"PseudoHaz <- predict(Qj.mod, type='response', newdata=x)",sep="")))
          if(j != min(J)){
            eval(parse(text=paste("x$hazLessThan",j," <- rowSums(cbind(rep(0, nrow(x)),x[,paste0('Q',J[J<j],'Haz')]))",sep="")))
            eval(parse(text=paste("x$Q",j,"Haz <- x$Q",j,"PseudoHaz * (1-x$hazLessThan",j,")",sep="")))
          }else{
            eval(parse(text=paste("x$Q",j,"Haz <- x$Q",j,"PseudoHaz",sep="")))
          }
          x
        }, j=j)
      }
    }else{
      for(j in J){
        Qj.form <- sprintf("%s ~ %s", paste("N",j,sep=""), glmFormula.event)
        X <- model.matrix(as.formula(Qj.form),data=dataList[[1]])

        NlessthanJ <- rep(0, nrow(dataList[[1]]))
        for(i in J[J<j]){
          eval(parse(text=paste("NlessthanJ <- NlessthanJ + dataList[[1]]$N",i,sep="")))
        }

        dataList <- lapply(dataList, function(x,j){
          if(j != min(J)){
            eval(parse(text=paste("x$hazLessThan",j," <- rowSums(cbind(rep(0, nrow(x)),x[,paste0('Q',J[J<j],'Haz')]))",sep="")))
          }else{
            eval(parse(text=paste("x$hazLessThan",j," <- 0",sep="")))
          }
          x
        },j=j)

        eval(parse(text=paste("Ytilde <- (dataList[[1]]$N",j,"-dataList[[1]]$l",j,")/(pmin(dataList[[1]]$u",j,", 1 - dataList[[1]]$hazLessThan",j,")  - dataList[[1]]$l",j,")",sep="")))
        fm <- optim(par=rep(0,ncol(X)), fn=LogLikelihood, Y=Ytilde, X=X,
                    method="BFGS",gr=grad,
                    control=list(reltol=tol/1e4))
        if(fm$convergence!=0){
          return("convergence failure")
        }else{
          beta <- fm$par

          dataList <- lapply(dataList, function(x,j){
            newX <- model.matrix(as.formula(Qj.form),data=x)
            eval(parse(text=paste("x$Q",j,"PseudoHaz <- plogis(newX%*%beta)",sep="")))
            eval(parse(text=paste("x$Q",j,"Haz <- (pmin(x$u",j,", 1 - x$hazLessThan",j,")  - x$l",j,") * x$Q",j,"PseudoHaz + x$l",j,sep="")))
            x
          },j=j)
        }
      }
    }
  }else if(is.null(glmFormula.event) & is.null(strata)){
    for(j in J){
      # add up all events less than current j to see who to include in regression
      NlessthanJ <- rep(0, nrow(dataList[[1]]))
      for(i in J[J<j]){
        eval(parse(text=paste("NlessthanJ <- NlessthanJ + dataList[[1]]$N",i,sep="")))
      }

      Qj.mod <- eval(parse(text=paste("SuperLearner(Y=dataList[[1]]$N",j,"[NlessthanJ==0],
                             X=dataList[[1]][NlessthanJ==0,c('t', 'Z', names(adjustVars))],
                             id=dataList[[1]]$id[NlessthanJ==0],
                             family=binomial(),
                             SL.library=SLlibrary.event,
                             verbose=verbose)",sep="")))

      # get predictions back
      dataList <- lapply(dataList, function(x,j){
        eval(parse(text=paste("x$Q",j,"PseudoHaz <- predict(Qj.mod, onlySL=T, newdata=x)[[1]]",sep="")))
        if(j != min(J)){
          eval(parse(text=paste("x$hazLessThan",j," <- rowSums(cbind(rep(0, nrow(x)),x[,paste0('Q',J[J<j],'Haz')]))",sep="")))
          eval(parse(text=paste("x$Q",j,"Haz <- x$Q",j,"PseudoHaz * (1-x$hazLessThan",j,")",sep="")))
        }else{
          eval(parse(text=paste("x$Q",j,"Haz <- x$Q",j,"PseudoHaz",sep="")))
        }
        x
      }, j=j)
    }
  }else if(is.null(glmFormula.event) & is.null(SLlibrary.event)){
    for(j in J){
      X <- model.matrix(~factor(dataList[[1]]$t)*dataList[[1]]$Z*factor(dataList[[1]]$strata))
      eval(parse(text=paste("beta <- solve(t(X)%*%X)%*%t(X)%*%dataList[[1]]$N",j,sep="")))
      dataList <- lapply(dataList, function(x,beta){
        X <- model.matrix(~factor(x$t)*x$Z*factor(x$strata))
        eval(parse(text=paste("x$Q",j,"Haz <- X%*%beta",sep="")))
        x
      },beta=beta)
    }
  }
  dataList

}




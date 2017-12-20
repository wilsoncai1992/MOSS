#' Validate and preprocess the data
#'
#' @param      dat   The dat
#' @param      dW    The d w
#' @param      nbin  how many levels should continuous W's be binned into
#'
#' @return
#' @export
#'
check_and_preprocess_data <- function(dat, dW, nbin = 4, T.cutoff = NULL) {
  to_keep <- (dat$T.tilde != 0)
  dW <- dW[to_keep]
  dat <- dat[to_keep,]
  
  n.data <- nrow(dat)
  W_names <- grep('W', colnames(dat), value = TRUE)
  # ==================================================================================
  # input validation
  # ==================================================================================
  if (length(dW) != n.data) stop('The length of input dW is not same as the sample size!')
  
  if (!('delta' %in% colnames(dat))) {
    warning('delta not found. Set delta = 1.')
    dat$delta <- rep(1, nrow(dat))
  }
  
  if ('T.TILDE' %in% toupper(colnames(dat))) {
    # if there is t.tilde in variable, remove any T
    keeps <- setdiff(colnames(dat), 'T')
    dat <- dat[,keeps]
  }else if('T' %in% toupper(colnames(dat))){
    message('no t.tilde, rename T to T.tilde')
    dat <- rename(dat, T.tilde = T)
  }else{
    # if there are no T.tilde
    stop("There should be T.tilde variable!")
  }
  # ==================================================================================
  # continuous W: binning
  # ==================================================================================
  checkBinary <- function(v, naVal="NA") {
    if (!is.numeric(v)) stop("Only numeric vectors are accepted.")
    
    vSet = unique(v)
    if (!missing(naVal)) vSet[vSet == naVal] = NA
    vSet = vSet[!is.na(vSet)]
    
    if (any(as.integer(vSet) != vSet)) "con"
    else if (length(vSet) > 2) "con"
    else "bin"
  }
  checkBinary_df <- function(df, W_names) {
    all_out <- character()
    for (it in W_names) {
      all_out <- c(all_out, checkBinary(v = df[,it]))
    }
    return(all_out)
  }
  is_conti <- checkBinary_df(df = dat, W_names = W_names) != 'bin'
  W_conti <- W_names[is_conti]
  
  if(nbin <=0) {
    message('Not binning any variable')
  }else{
    if(any(is_conti)) message(paste("Binning continuous covariates:", W_conti, collapse = ','))
    for (it in W_conti) {
      dat[,it] <- as.numeric(cut(dat[,it], nbin))
    }
  }
  # =======================================================================================
  # user-applied cutoff
  # =======================================================================================
  if(!is.null(T.cutoff)){
    message(paste('Manual right-censor the data up to time', T.cutoff))
    dat$delta[dat$T.tilde >= T.cutoff] <- 0
    dat$T.tilde[dat$T.tilde >= T.cutoff] <- T.cutoff
  }
  # ==================================================================================
  # check for positivity
  # ==================================================================================
  check_positivity(dat = dat, posit_level = 0.05)
  
  return(list(dat = dat, dW = dW, n.data = n.data, W_names = W_names))
}


#' Check positivity assumption of input data
#'
#' Check whether positivity assumption is violated in data
#' censored data NOT supported
#'
#' @param dat input data.frame
#' @param posit_level level of positivity to detect from data
#'
#' @return raise warning when positivity violated
#' @export
#'
#' @examples
check_positivity <- function(dat, posit_level = 0.05) {
  W_names <- grep('W', colnames(dat), value = TRUE)
  W <- dat[,W_names]
  A <- dat[,'A']
  
  table_out <- ftable(cbind(W,A))
  # if W is vector
  if (is.null(dim(W))) table_out <- table(data.frame(cbind( W,A)))
  table_out <- table_out/rowSums(table_out)
  # there are strata where there is NO observation
  table_out <- table_out[complete.cases(table_out),]
  # table_out[is.na(table_out)] <- 0
  
  posit_violate_by_strata <- rowSums(table_out < posit_level)
  # some strata have absolutely no entry
  # posit_violate_by_strata[posit_violate_by_strata == 2] <- 1
  
  if(any(as.logical(posit_violate_by_strata))) warning(paste('Positivity assumption violated', sum(posit_violate_by_strata), 'out of', length(posit_violate_by_strata), 'stratas!'))
}

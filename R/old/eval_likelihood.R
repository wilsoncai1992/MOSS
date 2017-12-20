# #' Title
# #'
# #' @param obj
# #'
# #' @return
# #' @export
# #'
# #' @examples
# #' #NA
# eval_loglike <- function(obj) {
#     Qn.current_full <- obj$Qn.current_full
#     h.hat.t_full_current <- obj$h.hat.t_full_current

#     delta <- obj$dat$delta
#     T.tilde <- obj$dat$T.tilde

#     # log components
#     logS <- log(Qn.current_full)
#     logh <- numeric()
#     for(it in 1:nrow(obj$dat)){
#         if((T.tilde[it]+1) <= ncol(Qn.current_full)){
#             logS[it,(T.tilde[it]+1) : ncol(Qn.current_full)] <- 0
#         }
#         logh[it] <- log(h.hat.t_full_current[it,T.tilde[it]])
#     }

#     # log-likelihood for each n
#     loglike <- rowSums(logS) + delta * logh
#     # TEMP: remove those not coorespond to counterfactual survival
#     # loglike <- loglike[obj$dat$A == 0]
#     # loglike <- loglike[obj$dat$A == 1]

#     loglike <- sum(loglike)

#     return(loglike)
# }


#' Title
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
#' #NA
eval_loglike <- function(obj, dW) {
  
  Qn.current_full <- obj$Qn.current_full
  h.hat.t_full_current <- obj$h.hat.t_full_current
  
  delta <- obj$dat$delta
  T.tilde <- obj$dat$T.tilde
  
  # log components
  logS <- numeric()
  logh <- numeric()
  for(it in 1:nrow(obj$dat)){
    logS[it] <- log(Qn.current_full[it,T.tilde[it]])
    logh[it] <- log(h.hat.t_full_current[it,T.tilde[it]])
  }
  
  # log-likelihood for each n
  loglike <- logS + delta * logh
  # TEMP: remove those not coorespond to counterfactual survival
  loglike <- loglike[obj$dat$A == dW]
  # loglike <- loglike[obj$dat$A == 1]
  
  loglike <- sum(loglike)
  
  return(loglike)
}


#' Title
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
#' #NA
plot_loglike_path <- function(obj) {
  all_loglikeli <- obj$convergence$all_loglikeli
  all_stopping <- obj$convergence$all_stopping
  par(mar=c(5,4,4,5)+.1)
  plot(all_loglikeli, type = 'l', ylab = 'log-likelihood', xlab = 'epoch')
  par(new=TRUE)
  plot(all_stopping,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
  axis(4)
  mtext("L2 norm of P_n(EIC)",side=4,line=3)
  legend("right",col=c("black","blue"),lty=1,legend=c("log-likelihood","L2 norm of P_n(EIC)"))
}

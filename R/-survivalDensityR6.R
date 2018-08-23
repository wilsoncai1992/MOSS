# #' @export
# survivalDensity <- R6Class("survivalDensity",
#   public = list(
#     pdf = NULL,
#     #
#     S = NULL,
#     h = NULL,
#     initialize = function(pdf){
#       if(class(pdf) == 'discreteDensity') {
#         self$pdf <- pdf
#       }else{
#         stop('type error')
#       }
#     },
#     compute_survival_from_pdf = function(){
#       self$S <- compute_step_cdf(pdf.mat = self$pdf$p, t.vec = 1:max(self$pdf$t_grid), start = Inf)
#     },
#     compute_hazard_from_pdf_and_survival = function(){
#       hazard_new <- matrix(0, nrow = nrow(self$pdf), ncol = max(self$pdf$t_grid))
#       for (it in 1:nrow(self$pdf)) {
#         hazard_new[it, ] <- self$pdf$p[it, ] / self$S[it,]
#       }
#       # dirty fix: upper bound hazard
#       hazard_new[hazard_new >= 1] <- .8

#       self$h <- hazard_new
#     },
#     plot_survival_curve = function(...){
#       step_curve <- stepfun(x = 1:max(self$pdf$t_grid), y = c(1, colMeans(self$S)))
#       # can `add`, `col`
#       curve(step_curve, from = 0, to = max(self$pdf$t_grid), ...)
#     }
#   )
# )

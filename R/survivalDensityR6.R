#' @export
survivalDensity <- R6Class("survivalDensity",
  public = list(
    pdf = NULL,
    #
    P = NULL,
    initialize = function(pdf){
      self$pdf <- pdf
    },
    compute_survival_from_pdf = function(){
      compute_step_cdf(pdf.mat = self$pdf, t.vec = 1:max(self$pdf$t_grid), start = Inf)
    }
  )
)

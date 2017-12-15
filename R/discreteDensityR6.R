#' @export
discreteDensity <- R6Class("discreteDensity",
  public = list(
    p = NULL,
    t_grid = NULL,
    n = NULL,
    initialize = function(p, t_grid){
      self$p <- p
      self$t_grid <- t_grid
      self$n <- nrow(p)
    },
    normalize = function(true_integral = 1){
      self$p <- self$p / rowSums(self$p) * true_integral
    },
    multiply_vec_to_each_row = function(vec){
      self$p <- multiple_vector_to_matrix(self$p, vec)
    }
  )
)

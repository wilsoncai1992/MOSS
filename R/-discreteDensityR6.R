# #' @export
# discreteDensity <- R6Class("discreteDensity",
#   public = list(
#     p = NULL,
#     t_grid = NULL,
#     n = NULL,
#     initialize = function(p, t_grid = NULL){
#       self$p <- p
#       if (is.null(t_grid)) t_grid <- 1:ncol(p)
#       self$t_grid <- t_grid
#       self$n <- nrow(p)
#     },
#     normalize = function(true_integral = 1){
#       self$p <- self$p / rowSums(self$p) * true_integral
#     },
#     multiply_vec_to_each_row = function(vec){
#       self$p <- multiply_vector_to_matrix(self$p, vec)
#     },
#     plot = function(...){
#       step_curve <- stepfun(x = 1:max(self$t_grid), y = c(0, self$p))
#       # can `add`, `col`
#       curve(step_curve, from = 0, to = max(self$t_grid), ...)
#     }
#   )
# )

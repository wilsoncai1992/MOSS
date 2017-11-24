#' generate step function object using y value and corresponding jump points
#'
#' @param y.vec a vector of step function values
#' @param t.vec a vector for all jump locations of step function.
#' NOTE: the first element value of y.vec is the flat part LEFT of the first jump point specified in t.vec
#'
#' @return sfun a function object, that can perform any mapping
#' @export
#'
#' @examples
#' # TO DO
create_step_func <- function(y.vec, t.vec) {
  if (length(y.vec) != (length(t.vec) + 1)) warning('the legnth of input vectors incorrect!')
  sfun  <- stepfun(t.vec, y.vec, f = 0) # before the first jump point, the step function is 0 value
  return(sfun)
}

#' compute \eqn{I\{T.tilde >= t\}}
#'
#' loop over t.vec
#'
#' @param Time length n vector of failure time
#' @param t.vec t value of interest
#'
#' @return a binary vector, of length = t.vec
#' @export
#'
#' @examples
#' # TO DO
create_Yt_vector <- function(Time, t.vec) {
  (Time >= t.vec) + 0
}

#' compute \eqn{I\{T.tilde == t, Delta = 1\}}
#'
#' loop over t.vec
#'
#' @param Time length n vector of failure time
#' @param Delta length n vector of censoring indicator
#' @param t.vec t value of interest
#'
#' @return a binary vector, of length = t.vec
#' @export
#'
#' @examples
#' # TO DO
create_Yt_vector_with_censor <- function(Time, Delta, t.vec) {
  ((Time == t.vec) & (Delta == 1)) + 0
}

#' compute cumulative distribution function of a step-shaped (empirical) density
#'
#' @param pdf.mat if input vector = compute cdf for a step-function pdf;
#'              if input matrix = compute cdf for several step-function pdf with same jump points
#' @param t.vec unique jump points of step function
#' @param start if -Inf = from left to right; if Inf = from right to left.
#'
#' @return vector of cdf value
#' @export
#'
#' @examples
#' # TO DO
compute_step_cdf <- function(pdf.mat, t.vec, start = -Inf) {
  interval.size <- diff(t.vec)
  # interval.size <- c(0, interval.size)
  interval.size <- c(interval.size, 0) # 09-07

  # compute the mass
  if(is.matrix(pdf.mat)){
    # if input with multi-sample
    mass.by.interval <- sweep(pdf.mat,MARGIN=2, interval.size, `*`)
    # multiplies the interval length to each row of the y-values
    # the result is a matrix, each row is a single pdf, and entries are the mass

  }else{
    # if input with one-sample
    mass.by.interval <- pdf.mat * interval.size
  }


  if(is.infinite(start) & (start < 0)){
    # ======================================================================
    # start from -Inf
    if(is.matrix(pdf.mat)){
      # if input with multi-sample
      cdf.by.interval <- t(apply(mass.by.interval, 1, cumsum)) # cumsum of mass for each row, from left to right
    }else{
      # if input with one-sample
      cdf.by.interval <- cumsum(mass.by.interval)
    }
  }else{
    # ======================================================================
    # start from +Inf
    if(is.matrix(pdf.mat)){
      # if input with multi-sample
      cdf.by.interval <- t(apply(mass.by.interval, 1, function(obj) rev(cumsum(rev(obj))) ) )
    }else{
      # if input with one-sample
      cdf.by.interval <- rev(cumsum(rev(mass.by.interval)))
    }
  }
  # ======================================================================
  return(cdf.by.interval)
}


#' compute l2 inner product of two step functions
#'
#' f and g
#'
#' @param f.step two step-function pdf with shared jump points; can be matrix input: nrow = # of different step-function pdf, ncol = length(T.grid)
#' @param g.step two step-function pdf with shared jump points; can be matrix input: nrow = # of different step-function pdf, ncol = length(T.grid)
#' @param T.grid shared jump points
#'
#' @return scalar
#' @export
#'
#' @examples
#' # TO DO
l2_inner_prod_step <- function(f.step, g.step, T.grid) {
  if(is.vector(f.step) & is.vector(g.step)){
    # both f and g are one sample
    f.times.g <- f.step * g.step
  }
  if(!is.vector(f.step) & is.vector(g.step)){
    # f: multi-sample
    # g: one-sample
    f.times.g <- sweep(f.step,MARGIN=2,g.step,`*`) # multiply g to each row of f.
  }
  if(is.vector(f.step) & !is.vector(g.step)){
    # f: one-sample
    # g: multi-sample
    f.times.g <- sweep(g.step,MARGIN=2,f.step,`*`) # multiply f to each row of g.
  }
  if(!is.vector(f.step) & !is.vector(g.step)){
    # both f and g are multi-sample of same sample size
    if(nrow(f.step) != nrow(g.step)) stop('f and g have different sample size!')
    f.times.g <- f.step * g.step
  }
  # ------------------------------------------------------------------------------------
  result <- compute_step_cdf(f.times.g, T.grid)
  if(!is.vector(f.step) | !is.vector(g.step)){
    # there is multi-sample
    result <- apply(result, 1, function(obj) tail(obj, 1))
  }else{
    # both f and g are one-sample
    result <- tail(result, 1)
  }
  return(result)
}

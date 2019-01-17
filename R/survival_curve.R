library("R6")

#' @export
survival_curve <- R6Class("survival_curve",
  public = list(
    t = NULL,
    hazard = NULL,
    # P( T >= t)
    survival = NULL,
    pdf = NULL,
    initialize = function(t, hazard = NULL, survival = NULL, pdf = NULL) {
      # only supports integer grid
      from_hazard <- !is.null(hazard)
      from_survival <- !is.null(survival)
      from_pdf <- !is.null(pdf)
      if (from_hazard + from_survival + from_pdf > 1) {
        stop("cannot construct from both")
      }
      if (!all.equal(t, seq(range(t)[1], range(t)[2]))) {
        stop("t is not integer without gap")
      }
      self$t <- t
      if (from_hazard) {
        message("construct from hazard")
        if ("data.frame" %in% class(hazard)) hazard <- as.matrix(hazard)
        if ("numeric" %in% class(hazard)) hazard <- matrix(hazard, nrow = 1)
        self$hazard <- hazard
      }
      if (from_survival) {
        message("construct from survival")
        if ("data.frame" %in% class(survival)) survival <- as.matrix(survival)
        if ("numeric" %in% class(survival)) survival <- matrix(survival, nrow = 1)
        self$survival <- survival
      }
      if (from_pdf) {
        message("construct from pdf")
        if ("data.frame" %in% class(pdf)) pdf <- as.matrix(pdf)
        if ("numeric" %in% class(pdf)) pdf <- matrix(pdf, nrow = 1)
        self$pdf <- pdf
      }
    },
    n = function() {
      n1 <- nrow(self$hazard)
      n2 <- nrow(self$survival)

      return(ifelse(is.null(n1), n2, n1))
    },
    hazard_to_survival = function() {
      self$survival <- matrix(NA, nrow = self$n(), ncol = max(self$t))
      # for (i in 1:self$n()) {
      #   self$survival[i, ] <- cumprod(1 - self$hazard[i, ])
      # }

      # cumulative hazard approach
      # integral from t = 0 and not include right bound. e.g. [0, 1)
      hazard_integral <- cbind(0, self$hazard)
      hazard_integral <- hazard_integral[, -ncol(hazard_integral)]
      for (i in 1:self$n()) {
        self$survival[i, ] <- exp(- cumsum(hazard_integral[i, ]))
      }
      return(self)
    },
    hazard_to_pdf = function() {
      self$hazard_to_survival()
      # not good using the theory formula
      # self$pdf <- self$hazard * self$survival
      self$pdf <- matrix(NA, nrow = self$n(), ncol = max(self$t))
      for (i in 1:self$n()) {
        self$pdf[i, ] <- c(- diff(self$survival[i, ]), 0)
      }
      return(self)
    },
    pdf_to_survival = function() {
      pdf2 <- cbind(0, self$pdf)
      pdf2 <- pdf2[, -ncol(pdf2)]
      # transpose: so that each row is one curve
      self$survival <- 1 - t(apply(pdf2, 1, cumsum))
    },
    pdf_to_hazard = function() {
      self$pdf_to_survival()
      self$hazard <- self$pdf / self$survival
    },
    display = function(type, W = NULL) {
      library("ggplot2")
      if (is.null(W)) {
        df <- data.frame(t = rep(self$t, self$n()))
      } else {
        if (class(W) != "numeric") stop("W only be univariate vector")
        if (length(W) != self$n()) stop("W length not correct")
        # the first Tmax rows are for the first subject
        df <- data.frame(
          t = rep(self$t, self$n()),
          W = rep(W, each = length(self$t))
        )
      }
      if (type == "survival") {
        df$s <- as.vector(t(self$survival))
        if (!is.null(W)) {
          gg <- ggplot(df, aes(x = t, y = round(W, digits = 1), z = s)) +
            geom_raster(aes(fill = s), interpolate = TRUE) +
            xlim(c(1, max(self$t))) +
            ylab("W") +
            theme_bw()
        } else {
          gg <- ggplot(df, aes(x = t, y = s)) +
            geom_line() +
            theme_bw() +
            ylim(c(-.1, 1.1))
        }
      }
      if (type == "hazard") {
        df$hazard <- as.vector(t(self$hazard))
        if (!is.null(W)) {
          gg <- ggplot(df, aes(x = t, y = round(W, digits = 1), z = hazard)) +
            geom_raster(aes(fill = hazard), interpolate = TRUE) +
            xlim(c(1, max(self$t))) +
            ylab("W") +
            theme_bw()
        } else {
          gg <- ggplot(df, aes(x = t, y = hazard)) +
            geom_line() +
            theme_bw()
        }
      }
      if (type == "pdf") {}
      return(gg)
    }
  )
)

# S ~ A = 1, W, t plot

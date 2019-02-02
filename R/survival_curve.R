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
        # message("construct from hazard")
        if ("data.frame" %in% class(hazard)) hazard <- as.matrix(hazard)
        if ("numeric" %in% class(hazard)) hazard <- matrix(hazard, nrow = 1)
        self$hazard <- hazard
      }
      if (from_survival) {
        # message("construct from survival")
        if ("data.frame" %in% class(survival)) survival <- as.matrix(survival)
        if ("numeric" %in% class(survival)) survival <- matrix(survival, nrow = 1)
        self$survival <- survival
      }
      if (from_pdf) {
        # message("construct from pdf")
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
      # working
      self$survival <- matrix(NA, nrow = self$n(), ncol = max(self$t))
      for (i in 1:self$n()) {
        hazard_here <- c(0, self$hazard[i, ])
        hazard_here <- hazard_here[-length(hazard_here)]
        self$survival[i, ] <- cumprod(1 - hazard_here)
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
    },
    create_ggplot_df = function() {
      # only for marginal survival curve
      return(data.frame(t = self$t, s = as.numeric(self$survival)))
    },
    ci = function(
      A,
      T_tilde,
      Delta,
      density_failure,
      density_censor,
      g1W,
      psi_n,
      A_intervene,
      alpha = 0.05
    ) {
      eic_fit <- eic$new(
        A = A,
        T_tilde = T_tilde,
        Delta = Delta,
        density_failure = density_failure,
        density_censor = density_censor,
        g1W = g1W,
        psi = psi_n,
        A_intervene = A_intervene
      )$all_t(k_grid = self$t)
      sigma <- apply(eic_fit, 2, sd)
      lower <- psi_n - sigma * 1.96
      upper <- psi_n + sigma * 1.96
      return(data.frame(t = self$t, lower = lower, upper = upper))
    }
  )
)


#' @export
evaluate_metric <- R6Class("evaluate_metric",
  public = list(
    survival = NULL,
    survival_truth = NULL,
    initialize = function(survival = NULL, survival_truth = NULL) {
      # only work for a vector of survival probabilities
      self$survival <- survival
      self$survival_truth <- survival_truth
      return(self)
    },
    evaluate_cross_entropy = function() {
      l <- c()
      for (t in self$survival$t) {
        s <- self$survival$survival[t]
        s_truth <- self$survival_truth$survival[t]
        if (s_truth == 1) l[t] <- -log(s)
        if (s_truth == 0) l[t] <- -log(1 - s)
        if (s_truth > 0 & s_truth < 1) l[t] <- -(s_truth * log(s) + (1 - s_truth) * log(1 - s))
      }
      return(data.frame(t = self$survival$t, metric = l))
    },
    evaluate_mse = function() {
      l <- c()
      bias <- as.numeric(self$survival$survival - self$survival_truth$survival)
      mse <- as.numeric((self$survival$survival - self$survival_truth$survival) ^ 2)
      return(data.frame(
        t = self$survival$t,
        bias = bias,
        mse = mse,
        truth = self$survival_truth$survival
      ))
    }
  )
)

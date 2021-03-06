% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_curve.R
\docType{class}
\name{evaluate_metric}
\alias{evaluate_metric}
\title{evaluate one survival_curve against another as truth}
\format{
\code{\link{R6Class}} object.
}
\value{
Object of \code{\link{R6Class}} with methods
}
\description{
evaluate one survival_curve against another as truth

evaluate one survival_curve against another as truth
}
\section{Methods}{

evaluate_mse create data.frame that computes mse for each time point
}

\examples{
\donttest{
  evaluation <- evaluate_metric$new(survival, survival_truth)
  evaluation$evaluate_mse()
}
}
\keyword{data}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{survival}}{estimated survival}

\item{\code{survival_truth}}{true survival}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{evaluate_metric$new()}}
\item \href{#method-evaluate_cross_entropy}{\code{evaluate_metric$evaluate_cross_entropy()}}
\item \href{#method-evaluate_mse}{\code{evaluate_metric$evaluate_mse()}}
\item \href{#method-clone}{\code{evaluate_metric$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{evaluate_metric$new(survival = NULL, survival_truth = NULL)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-evaluate_cross_entropy"></a>}}
\if{latex}{\out{\hypertarget{method-evaluate_cross_entropy}{}}}
\subsection{Method \code{evaluate_cross_entropy()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{evaluate_metric$evaluate_cross_entropy()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-evaluate_mse"></a>}}
\if{latex}{\out{\hypertarget{method-evaluate_mse}{}}}
\subsection{Method \code{evaluate_mse()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{evaluate_metric$evaluate_mse()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{evaluate_metric$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

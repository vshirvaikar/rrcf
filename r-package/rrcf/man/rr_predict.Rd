% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/causal_forest.R
\name{rr_predict}
\alias{rr_predict}
\title{Predict with an RR causal forest}
\usage{
\method{rr_predict}(
  object,
  newdata = NULL,
  num.threads = NULL,
  ...
)
}
\arguments{
\item{object}{The trained forest.}

\item{newdata}{Points at which predictions should be made. Note
that this matrix should have the number of columns as the training
matrix, and that the columns must appear in the same order.}

\item{num.threads}{Number of threads used in training. If set to NULL, the software
automatically selects an appropriate amount.}

}
\value{
Vector of RELATIVE treatment effect predictions.
}
\description{
Gets estimates of RELATIVE tau(x) using a trained causal forest.
}

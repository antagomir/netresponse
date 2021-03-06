% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add.ellipse.R
\name{add.ellipse}
\alias{add.ellipse}
\title{Add ellipse to an existing plot}
\usage{
add.ellipse(
  centroid,
  covmat,
  confidence = 0.95,
  npoints = 100,
  col = "black",
  ...
)
}
\arguments{
\item{centroid}{Vector with two elements defining the ellipse centroid.}

\item{covmat}{Covariance matrix for the investigated data. Only diagonal
covariances supported.}

\item{confidence}{Confidence level determining the ellipse borders based on
the covariance matrix.}

\item{npoints}{Number of plotting points.}

\item{col}{Color.}

\item{...}{Other arguments to be passed.}
}
\value{
Used for plotting side effects.
}
\description{
Calculates and plots ellipse corresponding to specified
confidence interval in 2-dimensional plot
}
\author{
Leo Lahti \email{leo.lahti@iki.fi}
}
\keyword{utilities}

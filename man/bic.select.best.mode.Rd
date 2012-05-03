\name{bic.select.best.mode}
\alias{bic.select.best.mode}
\title{Description: Select optimal number of mixture components by adding components until
the increase in objective function is below threshold.}
\usage{
  bic.select.best.mode(x, max.modes, bic.threshold)
}
\arguments{
  \item{x}{dat vector (for univariate analysis) or a matrix
  (for multivariate analysis)}

  \item{max.modes}{Maximum number of modes to be checked
  for mixture model selection}

  \item{bic.threshold}{BIC threshold which needs to be
  exceeded before a new mode is added to the mixture.

  Returns:}
}
\value{
  Fitted latent class model (parameters and free energy)
}
\description{
  Arguments:
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse")
}
\keyword{utilities}


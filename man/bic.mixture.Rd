\name{bic.mixture}
\alias{bic.mixture}
\title{Description: Latent class analysis based on (infinite) Gaussian mixture model.
If the input is data matrix, a multivariate model is fitted; if the input is a vector, a univariate model is fitted}
\usage{
  bic.mixture(x, max.modes, bic.threshold = 0,
    min.modes = 1, ...)
}
\arguments{
  \item{x}{samples x features matrix for multivariate
  analysis, or a vector for univariate analysis}

  \item{max.modes}{Maximum number of modes to be checked
  for mixture model selection}

  \item{bic.threshold}{BIC threshold which needs to be
  exceeded before a new mode is added to the mixture.}

  \item{min.modes}{minimum number of modes}

  \item{...}{Further optional arguments to be passed

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


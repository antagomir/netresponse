\name{P.s.individual}
\alias{P.s.individual}
\title{Description: Probabiity density for individual sample
Mainly for internal use; documentation will be provided later. Tools for calculating densities with Gaussian mixture models.}
\usage{
  P.s.individual(dat, pars, log = TRUE)
}
\arguments{
  \item{dat}{features x samples data matrix for mixture
  modeling}

  \item{pars}{Gaussian mixture model parameters (diagonal
  covariances); list with elements mu (mean vectors), sd
  (covariance diagonals), w (weights). The mu and sd are
  component x features matrices, w is vector giving weight
  for each component.}

  \item{log}{Logical. Return densities in log domain.

  Returns:}
}
\value{
  Probability density
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
\keyword{internal}
\keyword{utilities}


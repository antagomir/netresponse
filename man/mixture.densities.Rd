\name{mixture.densities}
\alias{P.rS}       
\alias{P.r.s}       
\alias{P.rs.joint.individual}       
\alias{P.rs.joint}       
\alias{P.Sr}       
\alias{P.S}       
\alias{P.sr}       
\alias{P.s.individual}       

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Mixture model density estimation functions.}
\description{Functions for mixture model density estimation; mainly for
  internal use.}
\usage{
P.rS(dat, pars, log = TRUE)
P.r.s(dat, pars, log = TRUE)
P.rs.joint.individual(dat, pars, log = TRUE)
P.rs.joint(dat, pars, log = TRUE)
P.Sr(dat, pars, log = TRUE)
P.S(dat, pars, log = TRUE)
P.sr(dat, pars, log = TRUE)
P.s.individual(dat, pars, log = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{features x samples data matrix for mixture modeling}
  \item{pars}{Gaussian mixture model parameters (diagonal covariances); list with elements mu
    (mean vectors), sd (covariance diagonals), w (weights). The mu and
    sd are component x features matrices, w is vector giving weight
    for each component.}
  \item{log}{Logical. Return densities in log domain.}
}
\details{Mainly for internal use; documentation will be provided
  later. Tools for calculating densities with Gaussian mixture models.}
\value{Scalar, vector, or matrix of probability densities (or
  log-densities), depending on the function.}
\references{See citation("netresponse").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}
%\seealso{}
\examples{
# For internal use.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
\keyword{ utilities }

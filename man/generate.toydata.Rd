\name{generate.toydata}
\alias{generate.toydata}
\title{generate.toydata}
\description{Generate simulated data for testing purposes.}
\usage{
  generate.toydata(Dim = 3, Nc = 3, Ns = 200, sd0 = 3, rgam.shape = 2, rgam.scale = 2, rseed = 12346) 
}
\arguments{
  \item{Dim}{Dimensionality of data}
  \item{Nc}{Number of modes}
  \item{Ns}{Number of data points}
  \item{sd0}{Component spread}
  \item{rgam.shape}{Shape parameter for Gamma distribution}
  \item{rgam.scale}{Scale parameter for Gamma distribution}
}
\value{
  Simulated data matrix (samples x features), component means, and component sds. sample-to-component mapping.
}
\examples{
# li <- generate.toydata()
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse")
}
\keyword{internal}


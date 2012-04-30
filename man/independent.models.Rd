\name{independent.models}
\alias{independent.models}
\title{independent.models}
\usage{
  independent.models(datamatrix, params, mixture.method =
    "vdp")
}
\arguments{
  \item{datamatrix}{datamatrix}

  \item{params}{parameters}

  \item{mixture.method}{Specify the approach to use in
  mixture modeling. Options. vdp (nonparametric Variational
  Dirichlet process mixture model); bic (based on Gaussian
  mixture modeling with EM, using BIC to select the optimal
  number of components)}
}
\value{
  \item{nodes }{Model for each node} \item{C }{Costs for
  individual models}
}
\description{
  Mainly for internal use. Provide independent models for
  each node.
}
\examples{
#
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse")
}
\keyword{internal}


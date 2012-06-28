\name{remove.negative.edges}
\alias{remove.negative.edges}
\title{remove.negative.edges}
\usage{
  remove.negative.edges(network, delta, datamatrix, params)
}
\arguments{
  \item{network}{network}

  \item{delta}{associated cost function value changes for
  each node merge}

  \item{datamatrix}{datamatrix}

  \item{params}{parameters}
}
\value{
  Filtered network
}
\description{
  Remove edges where the connected nodes correlate
  negatively
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


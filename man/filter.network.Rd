\name{filter.netw}
\alias{filter.netw}
\title{filter.netw}
\description{Mostly for internal use. Prefilter edges if speedups required.}
\usage{
filter.netw(network, delta, datamatrix, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{network}
  \item{delta}{associated cost function value changes for each node merge}
  \item{datamatrix}{datamatrix}
  \item{params}{parameters}
}
\details{Include to the network only the edges with the highest mutual information, calculated based on the first principal components.}
\value{Filtered network}
\references{
See citation("netresponse")
}
\author{Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
        Maintainer: Leo Lahti \email{leo.lahti@iki.fi}}
\examples{
# TBA
}
\keyword{ internal }


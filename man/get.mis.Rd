\name{get.mis}
\alias{get.mis}
\title{get.mis}
\description{Mainly for internal use. Estimate mutual information for node pairs based on the first principal components}
\usage{
get.mis(datamatrix, network, delta, network.nodes, G, params)
}
\arguments{
  \item{datamatrix}{datamatrix}
  \item{network}{network}
  \item{delta}{delta}
  \item{network.nodes}{network.nodes}
  \item{G}{G}
  \item{params}{params}
}

\value{
mutual information matrix
}

\references{
See citation("netresponse")
}
\author{Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
        Maintainer: Leo Lahti \email{leo.lahti@iki.fi}}
\examples{
# TBA
}
\keyword{ internal }

\name{pick.model.pairs}
\alias{pick.model.pairs}
\title{pick.model.pairs}
\description{Mainly for internal use. Calculate joint models for each node pair}
\usage{
pick.model.pairs(network, network.nodes, model.nodes, datamatrix, params)
}
\arguments{
  \item{network}{network}
  \item{network.nodes}{network.nodes}
  \item{model.nodes}{model.nodes}
  \item{datamatrix}{datamatrix}
  \item{params}{parameters}
}

\value{
  \item{model.pairs}{joint models for each node pair}
  \item{delta}{corresponding delta value for the cost function}
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


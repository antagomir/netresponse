\name{pick.model.pairs}
\alias{pick.model.pairs}
\title{pick.model.pairs}
\usage{
  pick.model.pairs(network, network.nodes, model.nodes,
    datamatrix, params, mixture.method = "vdp")
}
\arguments{
  \item{network}{network}

  \item{network.nodes}{network.nodes}

  \item{model.nodes}{model.nodes}

  \item{datamatrix}{datamatrix}

  \item{params}{parameters}

  \item{mixture.method}{Specify the approach to use in
  mixture modeling. Options. vdp (nonparametric Variational
  Dirichlet process mixture model); bic (based on Gaussian
  mixture modeling with EM, using BIC to select the optimal
  number of components)}
}
\value{
  \item{model.pairs}{joint models for each node pair}
  \item{delta}{corresponding delta value for the cost
  function}
}
\description{
  Mainly for internal use. Calculate joint models for each
  node pair
}
\examples{
#
}
\author{
  Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
  Maintainer: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse")
}
\keyword{internal}


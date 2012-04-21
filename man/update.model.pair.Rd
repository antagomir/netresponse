\name{update.model.pair}
\alias{update.model.pair}
\title{update.model.pair}
\description{Mainly for internal use. Calculate joint model for given node pair and update delta accordingly.}
\usage{
update.model.pair(datamatrix, delta, network, edge, network.nodes, G, params, model.nodes, model.pairs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datamatrix}{datamatrix}
  \item{delta}{delta}
  \item{network}{network}
  \item{edge}{edge}
  \item{network.nodes}{network.nodes}
  \item{G}{G}
  \item{params}{params}
  \item{model.nodes}{model.nodes}
  \item{model.pairs}{model.pairs}
}

\value{
  \item{model.pairs }{model.pairs}
  \item{delta }{delta}
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

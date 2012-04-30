\name{plot.subnet}
\alias{plot.subnet}
\title{plot.subnet}
\usage{
  \method{plot}{subnet}(x, subnet.id, network, plot.names =
    TRUE, ...)
}
\arguments{
  \item{x}{Result from NetResponse (detect.responses
  function).}

  \item{subnet.id}{Subnet id.}

  \item{network}{Original network used in the modelling.}

  \item{plot.names}{Plot node names (TRUE) or indices
  (FALSE).}

  \item{...}{Further arguments for plot function.}
}
\value{
  Used for its side-effects. Returns a matrix that
  describes the investigated subnetwork.
}
\description{
  Plot the given subnetwork.
}
\examples{
#

# res <- detect.responses(D, netw, verbose = FALSE)
# net <- plot.subnet(res, subnet.idx = 1)
}
\author{
  Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
  Maintainer: Leo Lahti <leo.lahti@iki.fi>
}
\references{
  L. Lahti et al.: Global modeling of transcriptional
  responses in interaction networks. Submitted.
}
\keyword{utilities}


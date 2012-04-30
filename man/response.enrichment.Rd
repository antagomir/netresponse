\name{response.enrichment}
\alias{response.enrichment}
\title{Enrichment for a specified sample group in the given response.}
\usage{
  response.enrichment(subnet.id, model, s, response,
    method = "hypergeometric")
}
\arguments{
  \item{subnet.id}{Subnet.}

  \item{model}{NetResponseModel object.}

  \item{s}{User-defined sample group. For instance, samples
  belonging to a particular annotation class.}

  \item{response}{Response id (integer) within the subnet.}

  \item{method}{Enrichment method.}
}
\value{
  List with enrichment statistics, depending on enrichment
  method.
}
\description{
  Calculate enrichment values for a specified sample group
  in the given response.
}
\examples{
#
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse").
}
\seealso{
  order.responses
}
\keyword{utilities}


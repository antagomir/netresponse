\name{response.enrichment}
\alias{response.enrichment}
\title{Enrichment for a specified sample group in the given response.}
\usage{
  response.enrichment(qofz, annotation.sample, response,
    method = "hypergeometric")
}
\arguments{
  \item{qofz}{sample-mode assignment matrix}

  \item{annotation.sample}{User-defined sample group. For
  instance, samples belonging to a particular annotation
  class.}

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
#enr <- response.enrichment(subnet.id, models, sample, response, method)
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


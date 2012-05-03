\name{order.responses}
\alias{order.responses}
\title{order.responses}
\usage{
  order.responses(model, sample, method = "hypergeometric",
    min.size = 2, max.size = Inf, min.responses = 2,
    subnet.ids = NULL, verbose = FALSE)
}
\arguments{
  \item{model}{NetResponseModel object.}

  \item{sample}{Measure enrichment of this sample (set)
  across the observed responses.}

  \item{method}{'hypergeometric' measures enrichment of
  factor levels in this response; 'precision' measures
  response purity for each factor level; 'dependency'
  measures logarithm of the joint density between response
  and factor level vs. their marginal densities:
  log(P(r,s)/(P(r)P(s)))}

  \item{min.size,max.size,min.responses}{Optional
  parameters to filter the results based on subnet size and
  number of responses.}

  \item{subnet.ids}{Specify subnets for which the responses
  shall be ordered.  By default, use all subnets.}

  \item{verbose}{Follow progress by intermediate messages.}
}
\value{
  A data frame with elements 'ordered.responses' which
  gives a data frame of responses ordered by enrichment
  score for the investigated sample.  The subnetwork,
  response id and enrichment score are shown. The method
  field indicates the enrichment calculation method. The
  sample field lists the samples et for which the
  enrichments were calculated. The info field lists
  additional information on enrichment statistics.
}
\description{
  Orders the responses by association strength (enrichment
  score) to a given sample set. For instance, if the
  samples correspond to a particular experimental factor,
  this function can be used to prioritize the responses
  according to their association strength to this factor.
}
\note{
  Tools for analyzing end results of the model.
}
\examples{
#
# - for given sample/s (factor level), order responses (across all subnets) by association strength (enrichment score)
#order.responses(model, sample, method  = "hypergeometric") # overrepresentation
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse") for citation details.
}
\keyword{utilities}


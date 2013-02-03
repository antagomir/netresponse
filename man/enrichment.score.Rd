\name{enrichment.score}
\alias{enrichment.score}
\title{Enrichment for a specified sample group in the given response.}
\usage{
  enrichment.score(assignment.matrix, which.mode,
    annotation.samples, method = "hypergeometric",
    data = NULL)
}
\arguments{
  \item{assignment.matrix}{Samples-to-modes assignment
  matrix. Soft or hard.}

  \item{which.mode}{Mode for which the enrichment is
  calculated.}

  \item{annotation.samples}{Samples assigned to the class
  for which the enrichment is calculated}

  \item{method}{Enrichment method.}

  \item{data}{Optionally required by some methods: data
  (samples x features)}
}
\value{
  Enrichment score and info. Enrichment score varies from 0
  (low) to 1 (high enrichment).
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
\keyword{utilities}


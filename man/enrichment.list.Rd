\name{enrichment.list}
\alias{enrichment.list}
\title{Description: Quantify association between modes and continuous variable}
\usage{
  enrichment.list(models, annotation.vector)
}
\arguments{
  \item{models}{List of models. Each model should have a
  sample-cluster assignment matrix qofz.}

  \item{annotation.vector}{annotation vector with discrete
  factor levels, and named by the samples

  Returns:}
}
\value{
  List with each element corresponding to one annotation
  variable and listing the responses according to
  association strength
}
\description{
  Arguments:
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse")
}
\keyword{utilities}


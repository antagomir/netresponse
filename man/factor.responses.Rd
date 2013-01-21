\name{factor.responses}
\alias{factor.responses}
\title{Description: List responses for each level of the given factor}
\usage{
  factor.responses(annotation.vector, model,
    method = "hypergeometric", min.size = 2, data = NULL)
}
\arguments{
  \item{annotation.vector}{annotation vector with discrete
  factor levels, and named by the samples}

  \item{model}{NetResponse model object}

  \item{method}{method for enrichment calculation}

  \item{min.size}{minimum sample size for a response}

  \item{data}{data (samples x features; or a vector in
  univariate case)

  Returns:}
}
\value{
  List with each element corresponding to one factor level
  and listing the responses according to association
  strength
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


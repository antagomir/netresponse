\name{list.responses}
\alias{list.responses}
\title{Description: List responses for all factors and levels in the given
annotation matrix}
\usage{
  list.responses(annotation.df, model,
    method = "hypergeometric", min.size = 2, qth = Inf,
    verbose = TRUE)
}
\arguments{
  \item{annotation.df}{annotation data.frame with discrete
  factor levels, rows named by the samples}

  \item{model}{NetResponse model object}

  \item{method}{method for enrichment calculation}

  \item{min.size}{minimum sample size for a response}

  \item{qth}{q-value threshold}

  \item{verbose}{verbose Returns:}
}
\value{
  Table listing all associations between the factor levels
  and responses
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


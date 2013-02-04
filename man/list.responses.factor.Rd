\name{list.responses.factor}
\alias{list.responses.factor}
\title{Description: List responses for all factors and levels in the given
annotation matrix}
\usage{
  list.responses.factor(annotation.df, models,
    method = "hypergeometric", min.size = 2, qth = Inf,
    verbose = TRUE, data = NULL, rounding = NULL)
}
\arguments{
  \item{annotation.df}{annotation data.frame with discrete
  factor levels, rows named by the samples}

  \item{models}{List of models. Each model should have a
  sample-cluster assignment matrix qofz.}

  \item{method}{method for enrichment calculation}

  \item{min.size}{minimum sample size for a response}

  \item{qth}{q-value threshold}

  \item{verbose}{verbose}

  \item{data}{data (samples x features; or a vector in
  univariate case)}

  \item{rounding}{rounding digits

  Returns:}
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


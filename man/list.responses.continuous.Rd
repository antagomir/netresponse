\name{list.responses.continuous}
\alias{list.responses.continuous}
\title{Description: Investigate association of a continuous variable and the modes}
\usage{
  list.responses.continuous(annotation.df, groupings,
    method = "t-test", qth = Inf, verbose = TRUE,
    rounding = NULL)
}
\arguments{
  \item{annotation.df}{annotation data.frame with discrete
  factor levels, rows named by the samples}

  \item{groupings}{List of mode groupings. Each element
  lists the samples assignment matrix qofz, or a vector of
  cluster indices named by the samples.}

  \item{method}{method for quantifying the association}

  \item{qth}{q-value threshold}

  \item{verbose}{verbose}

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


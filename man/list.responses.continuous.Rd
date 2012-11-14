\name{list.responses.continuous}
\alias{list.responses.continuous}
\title{Description: Investigate association of a continuous variable and the modes}
\usage{
  list.responses.continuous(annotation.df, model,
    method = "t-test", min.size = 1, qth = Inf,
    verbose = TRUE, data = NULL)
}
\arguments{
  \item{annotation.df}{annotation data.frame with discrete
  factor levels, rows named by the samples}

  \item{model}{NetResponse model object}

  \item{method}{method for quantifying the association}

  \item{min.size}{minimum sample size for a response}

  \item{qth}{q-value threshold}

  \item{verbose}{verbose}

  \item{data}{data}

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


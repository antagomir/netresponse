\name{plot.data}
\alias{plot.data}
\title{Plot observed data.}
\usage{
  \method{plot}{data}(x, subnet.id, labels, ...)
}
\arguments{
  \item{x}{NetResponseModel object.}

  \item{subnet.id}{Specify the subnetwork.}

  \item{labels}{Annotation categories.}

  \item{...}{Further arguments for plot function.

  Return:}
}
\value{
  ggplot2 plot object
}
\description{
  Plotting tool for measurement data.
}
\details{
  Produces boxplot for each feature in each annotation
  category for the selected subnetwork.
}
\examples{
#
}
\author{
  Leo Lahti <leo.lahti@iki.fi>
}
\references{
  See citation("netresponse")
}
\seealso{
  plot.responses
}
\keyword{utilities}


\name{plot.data}
\alias{plot.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot observed data.}
\description{Plotting tool for measurement data.}
%\usage{plot.data(x, subnet.id, labels, ...)}
\usage{\method{plot}{data}(x, subnet.id, labels, \dots)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{NetResponseModel object.}
  \item{subnet.id}{Specify the subnetwork.}
  \item{labels}{Annotation categories.}
  \item{...}{Further arguments for plot function.}
}
\details{Produces boxplot for each feature in each annotation category for the selected subnetwork.}
\value{ggplot2 plot object}
\references{See citation("netresponse")}
\author{Leo Lahti <leo.lahti@iki.fi>}
\seealso{plot.responses}
\examples{
## plot.data(x, subnet.id, labels)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }


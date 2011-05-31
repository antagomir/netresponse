\name{plot.associations}
\alias{plot.associations}
\title{Association strength between category labels and responses.}
\description{Plot association strength between user-defined category labels and responses in a selected subnetwork.}
\usage{
plot.associations(x, subnet.id, labels, method = "hypergeometric", mode = "group.by.classes", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{NetResponseModel object}
  \item{subnet.id}{Subnetwork.}
  \item{labels}{Factor. Labels for the data samples. Name by samples, or provide in the same order as in the original data.}
  \item{method}{Method to calculate association strength.}
  \item{mode}{group.by.responses or group.by.classes: indicate barplot grouping type.}
  \item{\dots}{Other arguments to be passed for plot.}
}
\details{Associations are showm in terms -log10(p) enrichment values for the annotation categories for the responses within the specified subnetwork. No correction for multiple testing.}
\value{Used for side effect (plotting).}
\references{See citation("netresponse").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}
\seealso{plot.responses}
\examples{
#
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }


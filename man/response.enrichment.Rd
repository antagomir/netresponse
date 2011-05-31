\name{response.enrichment}
\alias{response.enrichment}
\title{Enrichment for a specified sample group in the given response.}
\description{Calculate enrichment values for a specified sample group in the given response.}
\usage{
response.enrichment(subnet.id, model, s, response, method = "hypergeometric")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subnet.id}{Subnet.}
  \item{model}{NetResponseModel object.}
  \item{s}{User-defined sample group. For instance, samples belonging to a particular annotation class.}
  \item{response}{Response id (integer) within the subnet.}
  \item{method}{Enrichment method.}
}
%\details{}
\value{List with enrichment statistics, depending on enrichment method.}
\references{See citation("netresponse").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{order.responses}
\examples{
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }


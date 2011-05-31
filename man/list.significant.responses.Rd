\name{list.significant.responses}
\alias{list.significant.responses}
\title{List responses with significant associations to a given sample group.}
\description{List responses with significant associations to a given sample group.}
\usage{
list.significant.responses(model, sample, qth = 1, method = "hypergeometric")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{NetResponseModel object.}
  \item{sample}{User-specified samples group for which the enrichments are calculated. For instance, an annotation category.}
  \item{qth}{q-value threshold for enrichments}
  \item{method}{Enrichment method.}
}
%\details{}
\value{Table containing statistics of the significantly associated responses.}
\references{See citation("netresponse")}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}
\seealso{response.enrichment}
\examples{
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }


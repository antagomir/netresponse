\name{order.responses}
\alias{order.responses}
\title{
 order.responses
}
\description{
 Orders the responses by association strength (enrichment score) to a
given factor level.
}
\usage{
%order.responses(model, phenodata, which.factor, which.level, method = "hypergeometric")
order.responses(model, sample, method = "hypergeometric")
}
\arguments{
  \item{model}{ NetResponseModel object. }
%  \item{phenodata}{ phenoData as in AffyBatch objects. A data frame of samples x factors.}
%  \item{which.factor}{ Specify the factor to investigate. }
%  \item{which.level}{ Which factor level to investigate. }
  \item{sample}{Measure enrichment of this sample (set) across the
observed responses.}
  \item{method}{ 'hypergeometric' measures enrichment of factor levels in this
response; 'precision' measures response purity for each factor level;
'dependency' measures logarithm of the joint density between response
and factor level vs. their marginal densities: log(P(r,s)/(P(r)P(s)))}
}
%\details{}
\value{
A data frame with elements 'ordered.responses' which gives a data
frame of responses ordered by enrichment score for the investigated
sample. The subnetwork, response id and enrichment score are
shown. The method field
indicates the enrichment calculation method. The sample field lists
the samples et
for which the enrichments were calculated.
}
\references{ See citation("netresponse") for citation details. }
\author{ Leo Lahti \email{leo.lahti@iki.fi}}
\note{ Tools for analyzing end results of the model. }
%\seealso{\code{\link{order.samples}}}
\examples{
# - for given sample/s (factor level), order responses (across all subnets) by association strength (enrichment score)
#order.responses(model, sample, method  = "hypergeometric") # overrepresentation

}
\keyword{ utilities }



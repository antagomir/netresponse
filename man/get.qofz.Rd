\name{getqofz}
\alias{getqofz}
\alias{getqofz,NetResponseModel-method}
\title{Sample-to-response matrix of probabilities P(r|s).}
\description{Retrieve P(r|s) from NetResponseModel model.}
\usage{
getqofz(model, subnet.id, log = FALSE)
}
\arguments{
  \item{model}{NetResponseModel object.}
  \item{subnet.id}{Subnetwork to investigate.}
  \item{log}{Output in log probabilities.}
}
\details{Calculates probability density for each response on a given
sample based on the estimated Gaussian mixture model.}
\value{Samples x responses matrix. Each entry is a probability P(r|s).}
\references{See citation("netresponse").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}
%% ~Make other sections like Warning with \section{Warning }{....} ~
%\seealso{}
\examples{
# qofz <- getqofz(model, subnet.id, log = FALSE)
}
\keyword{ utilities }
\keyword{ internal }

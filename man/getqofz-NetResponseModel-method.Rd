% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NetResponseModel-accessors.R
\name{getqofz,NetResponseModel-method}
\alias{getqofz,NetResponseModel-method}
\alias{getqofz}
\title{Sample-to-response matrix of probabilities P(r|s).}
\usage{
getqofz(model, subnet.id, log = FALSE)
}
\arguments{
\item{model}{NetResponseModel object.}

\item{subnet.id}{Subnetwork to investigate.}

\item{log}{Output in log probabilities.}
}
\value{
Samples x responses matrix. Each entry is a probability P(r|s).
}
\description{
Retrieve P(r|s) from NetResponseModel model.
}
\details{
Calculates probability density for each response on a given sample based on
the estimated Gaussian mixture model.
}
\examples{

# qofz <- getqofz(model, subnet.id, log = FALSE)

}
\references{
See citation('netresponse').
}
\author{
Leo Lahti \email{leo.lahti@iki.fi}
}
\keyword{internal}
\keyword{utilities}

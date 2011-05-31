\name{add.ellipse}
\alias{add.ellipse}
\title{Add ellipse to an existing plot.}
\description{Calculates and plots ellipse corresponding to specified confidence interval in 2-dimensional plot.}
\usage{
add.ellipse(centroid, covmat, confidence = 0.95, npoints = 100, col = "black", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{centroid}{Vector with two elements defining the ellipse centroid.}
  \item{covmat}{Covariance matrix for the investigated data. Only diagonal covariances supported.}
  \item{confidence}{Confidence level determining the ellipse borders based on the covariance matrix.}
  \item{npoints}{Number of plotting points.}
  \item{col}{Color.}
  \item{...}{Other arguments to be passed.}
}
%\details{}
\value{Used for plotting side effects.}
%\references{}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}
%\seealso{}
\examples{
#add.ellipse(centroid = c(0, 0), covmat = diag(c(1,2)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }


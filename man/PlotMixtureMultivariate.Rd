\name{PlotMixtureMultivariate}
\alias{PlotMixtureMultivariate}
\title{PlotMixtureMultivariate}
\usage{
  PlotMixtureMultivariate(x, means, sds, ws, labels = NULL, title = NULL, ...)
}
\arguments{
  \item{x}{data matrix (samples x features)}

  \item{means}{mode centroids (modes x features)}

  \item{sds}{mode standard deviations, assuming diagonal
  covariance matrices (modes x features, each row giving
  the sqrt of covariance diagonal for the corresponding
  mode)}

  \item{ws}{weight for each mode}

  \item{labels}{Optional: sample class labels to be
  indicated in colors.}

  \item{title}{title}

  \item{...}{Further arguments for plot function.}
}
\value{
  Used for its side-effects.
}
\description{
  Visualize data, centroids and response confidence
  intervals for a given Gaussian mixture model with PCA.
  Optionally, color the samples according to annotations
  labels.
}
\examples{
#plotMixture(dat, means, sds, ws)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse") for citation details.
}
\keyword{utilities}


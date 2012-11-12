\name{PlotMixtureUnivariate}
\alias{PlotMixtureUnivariate}
\title{PlotMixtureUnivariate}
\usage{
  PlotMixtureUnivariate(x, means, sds, ws, title.text = "", xlab.text = "", ylab.text = "", binwidth = 0.05, ...)
}
\arguments{

  \item{x}{data vector}

  \item{means}{mode centroids}

  \item{sds}{mode standard deviations}

  \item{ws}{weight for each mode}

  \item{title.text}{title text}

  \item{xlab.text}{x label text}

  \item{ylab.text}{y label text}

  \item{binwidth}{Histogram bin width}

  \item{...}{Further arguments to be passed}

}
\value{
  Used for its side-effects.
}
\description{
Visualize data, centroids and stds for a given univariate Gaussian mixture model.

}
\examples{
#PlotMixtureUnivariate(x, means, sds, ws)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse") for citation details.
}
\keyword{utilities}


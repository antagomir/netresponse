\name{PlotMixtureUnivariate}
\alias{PlotMixtureUnivariate}
\title{PlotMixtureUnivariate}
\usage{
  PlotMixtureUnivariate(x, means = NULL, sds = NULL,
    ws = NULL, title.text = NULL, xlab.text = NULL,
    ylab.text = NULL, binwidth = 0.05, qofz = NULL,
    density.color = "darkgray", cluster.assignments = NULL,
    ...)
}
\arguments{
  \item{x}{data vector}

  \item{means}{mode centroids}

  \item{sds}{mode standard deviations}

  \item{ws}{weight for each mode}

  \item{title.text}{Plot title}

  \item{xlab.text}{xlab.text}

  \item{ylab.text}{ylab.text}

  \item{binwidth}{binwidth for histogram}

  \item{qofz}{Mode assignment probabilities for each
  sample. Samples x modes.}

  \item{density.color}{Color for density lines}

  \item{cluster.assignments}{Vector of cluster indices,
  indicating cluster for each data point}

  \item{...}{Further arguments for plot function.

  Return:}
}
\value{
  Used for its side-effects
}
\description{
  Visualize data, centroids and stds for a given univariate
  Gaussian mixture model with PCA.
}
\details{
  Arguments:
}
\examples{
# plotMixtureUnivariate(dat, means, sds, ws)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse") for citation details.
}
\keyword{utilities}


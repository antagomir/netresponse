\name{plot.responses}
\alias{plot.responses}
\title{plot.responses}
\usage{
  \method{plot}{responses}(x, subnet.id, nc = 3, plot.names
    = TRUE, plot.mode = "network", xaxis = TRUE, yaxis =
    TRUE, plot.type = "twopi", mar = c(5, 4, 4, 2), horiz =
    TRUE, datamatrix = NULL, ...)
}
\arguments{
  \item{x}{Result from NetResponse (detect.responses
  function).}

  \item{subnet.id}{Subnet id.}

  \item{nc}{Number of columns for an array of images.}

  \item{plot.names}{Plot node names (TRUE) or indices
  (FALSE).}

  \item{plot.mode}{network: plot responses as a subnetwork
  graph; matrix, heatmap: plot subnetwork expression
  matrix. For both, expression of each gene is shown
  relative to the mean expression level of the gene;
  boxplot.data: feature-wise boxplots for hard
  sample-to-response assignments; response.barplot:
  estimated response centroids as barplot including 95%
  confidence intervals for the means; pca: PCA projection
  with estimated centroids and 95% intervals. In
  1-dimensional case a histogram is drawn. In
  two-dimensional case the original coordinates are used.}

  \item{xaxis,yaxis}{Logical. Plot row/column names.}

  \item{plot.type}{Network plot mode. For instance, 'neato'
  or 'twopi'.}

  \item{mar}{Figure margins.}

  \item{horiz}{Logical. Horizontal barplot.}

  \item{datamatrix}{Optionally provide the data matrix as argument.}

  \item{...}{Further arguments for plot function.}
}
\value{
  Used for its side-effects.
}
\description{
  Plot the detected transcriptional responses for a given
  subnetwork.
}
\examples{
#

#res <- detect.responses(D, netw)
#vis <- plot.responses(res, subnet.id)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse")
}
\seealso{
  \code{\link{plot.scale}}
}
\keyword{utilities}


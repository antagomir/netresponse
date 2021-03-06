% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot-methods.R
\name{plot_responses}
\alias{plot_responses}
\title{plot_responses}
\usage{
plot_responses(
  x,
  subnet.id,
  nc = 3,
  plot_names = TRUE,
  plot_mode = "network",
  xaxis = TRUE,
  yaxis = TRUE,
  plot_type = "twopi",
  mar = c(5, 4, 4, 2),
  horiz = TRUE,
  datamatrix = NULL,
  scale = FALSE,
  ...
)
}
\arguments{
\item{x}{Result from NetResponse (detect.responses function).}

\item{subnet.id}{Subnet id.}

\item{nc}{Number of columns for an array of images.}

\item{plot_names}{Plot node names (TRUE) or indices (FALSE).}

\item{plot_mode}{network: plot responses as a subnetwork graph; matrix,
heatmap: plot subnetwork expression matrix. For both, expression of each
gene is shown relative to the mean expression level of the gene;
boxplot_data: feature-wise boxplots for hard sample-to-response assignments;
response.barplot: estimated response centroids as barplot including 95%
confidence intervals for the means; pca: PCA projection with estimated
centroids and 95% intervals. In 1-dimensional case a histogram is drawn. In
two-dimensional case the original coordinates are used.}

\item{xaxis, yaxis}{Logical. Plot row/column names.}

\item{plot_type}{Network plot mode. For instance, 'neato' or 'twopi'.}

\item{mar}{Figure margins.}

\item{horiz}{Logical. Horizontal barplot_}

\item{datamatrix}{datamatrix}

\item{scale}{scale the phylotypes to unit length (only implemented for
plot_mode = 'matrix'}

\item{...}{Further arguments for plot function.}
}
\value{
Used for its side-effects.
}
\description{
Plot the detected transcriptional responses for a given
subnetwork.
plot_mode = 'network', xaxis = TRUE, yaxis = TRUE, plot_type = 'twopi', mar
= c(5, 4, 4, 2), horiz = TRUE, datamatrix = NULL, scale = FALSE, ...)
}
\examples{
#
#res <- detect.responses(D, netw)
#vis <- plot_responses(res, subnet.id)
}
\references{
See citation('netresponse')
}
\seealso{
\code{\link{plot_scale}}
}
\author{
Leo Lahti \email{leo.lahti@iki.fi}
}
\keyword{utilities}

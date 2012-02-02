\name{plot.scale}
\alias{plot.scale}
\title{plot.scale}
\description{ Plot the color scale used in visualization. }
\usage{
\method{plot}{scale}(x, y, m = NULL, cex.axis = 1.5, label.step = 2, interval = 0.1, two.sided = TRUE, label.start = 1, Nlab = 3, \dots)
}
\arguments{
  \item{x}{Breakpoints for the plot.}
  \item{y}{Color palette.}
  \item{m}{Breakpoints' upper limit.}
  \item{cex.axis}{Axis scale.}
  \item{label.step}{Density of the labels.}
  \item{interval}{Interval.}
  \item{two.sided}{Plot two-sided (TRUE) or one-sided (FALSE)
    visualization.}
  \item{label.start}{Label starting point.}
  \item{Nlab}{Number of labels to plot.}
  \item{...}{Further arguments for plot function.}
}
\value{Used for its side-effects.}
\references{See citation("netresponse")}
\author{Leo Lahti <leo.lahti@iki.fi>}
\note{Depends on Rgraphviz and igraph packages.}
\examples{
  #res <- detect.responses(D, netw, verbose = FALSE)
  #vis <- plot.responses(res, subnet.idx)
  #plot.scale(vis$breaks, vis$palette)
}
\keyword{utilities}

\name{plot.response}
\alias{plot.response}
\title{plot.response}
\description{Plot a specific transcriptional response for a given subnetwork.}
\usage{\method{plot}{response}(x, mynet, mybreaks, mypalette, plot.names = TRUE, colors = TRUE, plot.type = "twopi", \dots)}
%\usage{
%plot.response(
%x, mynet, mybreaks, mypalette, plot.names = TRUE, colors = TRUE, plot.type = "twopi", ...)
%}
\arguments{
  \item{x}{ A numerical vector, or NULL. }
  \item{mynet}{ Binary matrix specifying the interactions between
    nodes. }
  \item{mybreaks}{ Specify breakpoints for color plot. }
  \item{mypalette}{ Specify palette for color plot. }
  \item{plot.names}{ Plot node names (TRUE) or indices (FALSE).}
  \item{colors}{ Plot colors. Logical. }
  \item{plot.type}{ Network plot mode. For instance, 'neato' or 'twopi'.}
  \item{...}{Further arguments for plot function.}
}
\value{  Used for its side-effects. }
\references{L. Lahti et al.: Global modeling of transcriptional responses in interaction networks. Submitted.}
\author{Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
Maintainer: Leo Lahti <leo.lahti@iki.fi>}
\examples{
#tmp <- plot.response(model, mynet, maintext = paste("Subnetwork", subnet.id))
}
\keyword{utilities}

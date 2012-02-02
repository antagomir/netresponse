\name{plot.expression}
\alias{plot.expression}
\title{plot.expression}
\description{
Plot expression matrix in color scale. For one-channel data; plot
expression of each gene relative to its mean expression level over all
samples. Blue indicates decreased expression and red indicates
increased expression. Brightness of the color indicates magnitude of
the change. Black denotes no change.
}
%\usage{plot.expression(x, maintext, ...)}
\usage{\method{plot}{expression}(x, maintext, \dots)}
\arguments{
  \item{x}{samples x features matrix}
  \item{maintext}{main title}
  \item{\dots}{optional arguments}
}
%\details{}
\value{Used for its side effects.}
\references{See citation("netresponse").}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}
\seealso{\code{\link{plot.scale}}}
\examples{
#plot.expression(x)
}
\keyword{ utilities }


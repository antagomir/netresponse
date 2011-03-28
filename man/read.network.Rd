\name{read.network}
\alias{read.sif}
\title{Reading network files}
\description{Function to read network files.}
\usage{
read.sif(sif.file, format = "graphNEL", directed = FALSE)
}
\arguments{
  \item{sif.file}{Name of network file in SIF format.}
  \item{format}{Output format: igraph or graphNEL}
  \item{directed}{Logical. Directed/undirected graph. Not used in the current model.}
}
\details{Read in SIF network file, return R graph object in igraph or
  graphNEL format.}
\value{R graph object in igraph or
  graphNEL format.}
%\references{}
\author{Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}
%\seealso{}
\examples{
#net <- read.sif("network.sif")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }


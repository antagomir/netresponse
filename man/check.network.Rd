\name{check.network}
\alias{check.network}
\title{check.network}
\description{Mainly for internal use. Check input network and make formatting for detect.responses}
\usage{
check.network(network, datamatrix, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{Input network, see detect.responses}
  \item{datamatrix}{Input datamatrix, see detect.responses}
  \item{verbose}{Print intermediate messages}
}
\value{
  \item{formatted }{Formatted network (self-links removed)}
  \item{original }{Original network (possible in another representation format)}
  \item{delta }{Cost function changes corresponding to the 'formatted' network.}
  \item{nodes }{Nodes corresponding to the 'formatted' network.}
}
\references{
See citation("netresponse")
}
\author{Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
        Maintainer: Leo Lahti \email{leo.lahti@iki.fi}}
\seealso{detect.responses}
\examples{
# TBA
}
\keyword{ internal }


\name{get.subnets}
\alias{get.subnets}
\alias{get.subnets,NetResponseModel-method}
\title{get.subnets}
\arguments{
  \item{model}{Output from the detect.responses function.
  An object of NetResponseModel class.}

  \item{get.names}{Logical. Indicate whether to return
  subnetwork nodes using node names (TRUE) or node indices
  (FALSE).}

  \item{min.size,max.size}{Numeric. Filter out subnetworks
  whose size is not within the limits specified here.}

  \item{min.responses}{Numeric. Filter out subnetworks with
  less responses (mixture components) than specified here.}
}
\value{
  A list of subnetworks.
}
\description{
  List the detected subnetworks (each is a list of nodes in
  the corresponding subnetwork).
}
\examples{
# library(netresponse)
## Load a pre-calculated netresponse model obtained with
# model <- detect.responses(toydata$emat, toydata$netw, verbose = FALSE)
# data( toydata ); get.subnets(toydata$model)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  Leo Lahti et al.: Global modeling of transcriptional
  responses in interaction networks. Bioinformatics (2010).
  See citation("netresponse") for details.
}
\keyword{utilities}


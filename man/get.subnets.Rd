\name{get.subnets}
\alias{get.subnets}
\title{get.subnets}
\description{ List the detected subnetworks (each is a list of nodes
  in that subnetwork).
}
\usage{
subnets <- get.subnets(model, level = NULL, get.names = TRUE, stat = NULL, min.size = NULL, max.size = NULL, min.responses = NULL)
}
\arguments{

  \item{model}{Output from the detect.responses function. An object of
  NetResponseModel class.}

 \item{level}{ Agglomeration level to investigate. The agglomerative
    algorithm grows the subnetworks step-by-step. This option can be
    used to select a particular step during the learning process. Will
    be included in the next version. }

  \item{get.names}{Logical. Indicate whether to return subnetwork
  nodes using node names (TRUE) or node indices (FALSE).}

  \item{stat}{Subnetwork summary statistics. If this is not readily
  provided through this option (i.e. stat = NULL), it will be
  calculated. Can speed up the get.subnets function.}

  \item{min.size, max.size }{Numeric. Filter out subnetworks whose size
    is not within the limits specified here.}

  \item{min.responses }{Numeric. Filter out subnetworks with less
    responses (mixture components) than specified here.}

}
\value{ A list of subnetworks.}

\references{Leo Lahti et al.: Global modeling of transcriptional responses in interaction networks. Bioinformatics (2010).}

\author{Leo Lahti <leo.lahti@iki.fi>}

\examples{

library(netresponse)

# Load a pre-calculated netresponse model obtained with 
# model <- detect.responses(toydata$emat, toydata$netw, verbose = FALSE)
data( toydata )        
model <- toydata$model 

#List the detected subnetworks 
#(each is a list of nodes for the given subnetwork):
get.subnets(model)

}
\keyword{utilities}

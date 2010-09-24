\name{result.stats}
\alias{result.stats}
\title{result.stats}
\description{ Subnetwork statistics: size and number of distinct
  responses for each subnet. }
\usage{
result.stats( model, datamatrix, level )
}
\arguments{
  \item{model}{ Result from NetResponse (detect.responses function). }
  \item{datamatrix}{ Original data matrix for which the model was computed. }

  \item{level}{ Agglomeration level to investigate. The agglomerative
    algorithm grows the subnetworks step-by-step. This option can be
    used to select a specific step during the learning process. Will
    be included in the next version. }

}
\value{

  A 'subnetworks x properties' data frame containing the following
  elements.
  
  \item{subnet.size: }{ Vector of subnetwork sizes. }
  \item{subnet.responses: }{ Vector giving the number of responses in each subnetwork. }

}

\references{Leo Lahti et al.: Global modeling of transcriptional responses in interaction networks. Bioinformatics (2010).}

\author{Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
Maintainer: Leo Lahti <leo.lahti@iki.fi>}

\examples{

library(netresponse)

# Load a pre-calculated netresponse model obtained with 
# model <- detect.responses(toydata$emat, toydata$netw, verbose = FALSE)
data( toydata )        
model <- toydata$model # netresponse model
D <- toydata$emat         # data matrix

# Calculate summary statistics for the model
stat <- result.stats(model, D) 

}

\keyword{utilities}

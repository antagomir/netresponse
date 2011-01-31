\name{model.stats}
\alias{model.stats}
\title{model.stats}
\description{ Subnetwork statistics: size and number of distinct
  responses for each subnet. }
\usage{
model.stats( model )
}
\arguments{
  \item{model}{ Result from NetResponse (detect.responses function). }

%  \item{level}{ Agglomeration level to investigate. The agglomerative
%    algorithm grows the subnetworks step-by-step. This option can be
%    used to select a specific step during the learning process. Will
%    be included in the next version. }

}
\value{

  A 'subnetworks x properties' data frame containing the following
  elements.  
  \item{subnet.size: }{ Vector of subnetwork sizes. }
  \item{subnet.responses: }{ Vector giving the number of responses in each subnetwork. }

}

\references{Leo Lahti et al.: Global modeling of transcriptional
responses in interaction networks. Bioinformatics (2010).  See
citation("netresponse") for reference details.}

\author{Leo Lahti <leo.lahti@iki.fi>}

\examples{

library(netresponse)

# Load a pre-calculated netresponse model obtained with 
# model <- detect.responses(toydata$emat, toydata$netw, verbose = FALSE)
data( toydata )        
# Calculate summary statistics for the model
stat <- model.stats(toydata$model)

}

\keyword{utilities}

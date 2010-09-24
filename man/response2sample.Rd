\name{response2sample}
\alias{response2sample}
\title{response2sample}
\description{List the most strongly associated response of a given
  subnetwork for each sample.}
\usage{response2sample(model, subnet.id, datamatrix, component.list = TRUE)}
\arguments{

  \item{model}{A NetResponseModel object. Result from NetResponse
               (detect.responses function).}

  \item{subnet.id}{Subnet id. A natural number which specifies one of
  		   the subnetworks within the 'model' object.}

  \item{datamatrix}{Original data matrix used to obtain the 'model'
  		    object.}

  \item{component.list}{List samples separately for each mixture
                        component (TRUE). Else list the most strongly
                        associated component for each sample (FALSE).}

} 

\value{A list. Each element corresponds to one subnetwork response,
and contains a list of samples that are associated with the response
(samples for which this response has the highest probability
P(response | sample)).}

\references{Leo Lahti et al.: Global modeling of transcriptional responses in interaction networks. Bioinformatics (2010).}

\author{Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
Maintainer: Leo Lahti <leo.lahti@iki.fi>}
\examples{ 

library( netresponse )

# Load example data
data( toydata )         # Load toy data set
D    <- toydata$emat    # Response matrix (for example, gene expression)
model <- toydata$model  # Pre-calculated model

# Find the samples for each response (for a given subnetwork)
response2sample(model, subnet.id = 1, D)

}
\keyword{utilities}

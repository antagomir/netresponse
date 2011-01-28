\name{get.model.parameters}
\alias{get.model.parameters}
\title{get.model.parameters}
\description{Retrieve the mixture model parameters of the NetResponse
  algorithm for a given subnetwork.}
\usage{
get.model.parameters(model, subnet.id)
}
\arguments{
  \item{model}{Result from NetResponse (detect.responses function).}

 \item{ subnet.id }{Subnet identifier. A natural number which
  specifies one of the subnetworks within the 'model' object.}

% \item{level}{ Agglomeration level to investigate. The agglomerative
%    algorithm grows the subnetworks step-by-step. This option can be
%    used to select a specific step during the learning process. Will
%    be included in the next version. }

}
\value{
  A list with the following elements:
  \item{mu}{ Centroids for the mixture components. Components x nodes.}
  \item{sd}{ Standard deviations for the mixture components. A vector
    over the nodes for each component, implying the diagonal covariance
    matrix of the model (i.e. diag(std^2)). Components x nodes}
  \item{w}{Vector of component weights.}
  \item{nodes}{List of nodes in the subnetwork.}
  \item{K}{Number of mixture components.}
}
\details{Only the non-empty components are returned. Note: the original data matrix needs to be provided for function call separately.}

\references{Leo Lahti et al.: Global modeling of transcriptional
responses in interaction networks. Bioinformatics (2010). See
citation("netresponse") for details.}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\examples{
# Load toy data
data( toydata )          # Load toy data set
D     <- toydata$emat    # Response matrix (for example, gene expression)
model <- toydata$model   # Pre-calculated model

# Get model parameters for a given subnet
# (Gaussian mixture: mean, covariance diagonal, mixture proportions)
get.model.parameters(model, subnet.id = 1)

}

\keyword{utilities}

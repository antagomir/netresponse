\name{get.dat}
\alias{get.dat}
\title{get.dat}
\arguments{
  \item{model}{Result from NetResponse (detect.responses
  function).}

  \item{subnet.id}{Subnet identifier. A natural number
  which specifies one of the subnetworks within the 'model'
  object.}

  \item{sample}{Specify samples for which the data will be
  retrieved.}
}
\value{
  Data matrix features x samples.
}
\description{
  Retrieve data for a given subnetwork.
}
\examples{
# Load toy data
#data( toydata )          # Load toy data set
#D     <- toydata$emat    # Response matrix (for example, gene expression)
#model <- toydata$model   # Pre-calculated model
# Get model parameters for a given subnet
# (Gaussian mixture: mean, covariance diagonal, mixture proportions)
#get.dat(model, subnet.id = 1)
}
\author{
  Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse")
}
\keyword{utilities}


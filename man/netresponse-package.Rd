\name{netresponse-package}
\alias{netresponse-package}
\alias{netresponse}
\docType{package}
\title{NetResponse: Global modeling of transcriptional responses in interaction networks}
\description{Global modeling of transcriptional responses in interaction networks.}
\details{
\tabular{ll}{
Package: \tab netresponse\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2010-09-23\cr
License: \tab GNU GPL >=2\cr
LazyLoad: \tab yes\cr
}
}
\author{Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
Maintainer: Leo Lahti <leo.lahti@iki.fi>}
\references{Leo Lahti et al.: Global modeling of transcriptional responses in interaction networks. Bioinformatics (2010).}
\keyword{package}
%\seealso{%%~~ \code{\link[somepkg:somepkg-package]{somepkg}} ~~}}
\examples{
# Load the package
library(netresponse)

# Define parameters for toy data
Ns <- 200  # number of samples (conditions)
Nf <- 10   # number of features (nodes)
feature.names <- paste("feat", seq(Nf), sep="")
sample.names  <- paste("sample", seq(Ns), sep="") 

# random seed
set.seed( 123 )

# Random network
netw <- pmax(array(sign(rnorm(Nf^2)), dim = c(Nf, Nf)), 0)
# in pathway analysis nodes correspond to genes
rownames(netw) <- colnames(netw) <- feature.names

# Random responses of the nodes across conditions 
D <- array(rnorm(Ns*Nf), dim = c(Ns,Nf), dimnames = list(sample.names, feature.names))
D[1:100, 4:6]  <- t(sapply(1:(Ns/2),function(x){rnorm(3, mean = 1:3)}))
D[101:Ns, 4:6] <- t(sapply(1:(Ns/2),function(x){rnorm(3, mean = 7:9)}))

# Compute the model
model <- detect.responses(D, netw)

# Subnets (each is a list of nodes)
get.subnets( model )

# Retrieve model for the subnetwork with lowest cost function value
# means, standard devations and weights for the components
inds <- which(sapply(model@last.grouping, length) > 2)
subnet.id <- names(which.min(model@costs[inds]))
m <- get.model.parameters(model, subnet.id) 
print(m)
}

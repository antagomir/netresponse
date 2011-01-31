\name{sample2response}

\alias{sample2response}

\title{sample2response}

\description{ Probabilistic sample-response assignments for given
  subnet. }
\usage{sample2response(model, subnet.id)}
\arguments{
  \item{ model }{Result from NetResponse (detect.responses function).}

  \item{ subnet.id }{Subnet identifier. A natural number which
  specifies one of the subnetworks within the 'model' object.}

}

\value{A matrix of probabilities. Sample-response assignments for given
  subnet, listing the probability of each response, given a sample.}

\references{Leo Lahti et al.: Global modeling of transcriptional
responses in interaction networks. Bioinformatics (2010).  See
citation("netresponse") for citation details.}

\author{Leo Lahti \email{leo.lahti@iki.fi}}

\examples{
# Load toy data
data( toydata )        # Load toy data set
D    <- toydata$emat   # Response matrix (for example, gene expression)
netw <- toydata$netw   # Network

# Detect network responses
model <- detect.responses(D, netw, verbose = FALSE)

# Assign samples to responses (soft, probabilistic assignments sum to 1)
response.probabilities <- sample2response(model, subnet.id = 1)

}
\keyword{utilities}

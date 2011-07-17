\name{detect.responses}
\alias{detect.responses}
\title{detect.responses}
\description{ Main function of the NetResponse algorithm. Detecting network responses across the conditions. }
\usage{detect.responses(datamatrix, network, initial.responses = 1,
        max.responses = 10, max.subnet.size = 10, verbose =
        TRUE, prior.alpha = 1, prior.alphaKsi = 0.01, prior.betaKsi =
        0.01, update.hyperparams = 0, implicit.noise = 0, vdp.threshold =
        1.0e-5, merging.threshold = 0, ite = Inf,
        information.criterion = "AIC", speedup = TRUE,
        speedup.max.edges = 10, mc.cores = 1, ...)}

\value{NetResponseModel object.}

\references{Leo Lahti et al.: Global modeling of transcriptional
responses in interaction networks. Bioinformatics (2010). See
citation("netresponse") for details.}

\author{Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
        Maintainer: Leo Lahti \email{leo.lahti@iki.fi}}

\arguments{

  \item{datamatrix}{ Matrix of samples x features. For example, gene
    expression matrix with conditions on the rows, and genes on the
    columns. The matrix contains same features than the 'network'
    object, characterizing the network states across the different
    samples.

  }

  \item{network}{Network describing undirected pairwise interactions
    between features of 'datamatrix'. The following formats are
    supported: binary matrix, graphNEL, igraph, graphAM, Matrix,
    dgCMatrix, dgeMatrix}

  \item{initial.responses}{Initial number of components for each
  subnetwork model. Used to initialize calculations.}

  \item{max.responses }{Maximum number of responses for each
  subnetwork. Can be used to limit the potential number of network
  states.}

  \item{max.subnet.size}{Numeric. Maximum allowed subnetwork size.}

  \item{verbose}{Logical. Verbose parameter.}

  \item{implicit.noise}{Implicit noise parameter. Add implicit noise
  to vdp mixture model. Can help to avoid overfitting to local optima,
  if this appears to be a problem.}

  \item{update.hyperparams}{ Logical. Indicate whether to update
  hyperparameters during modeling.}

  \item{prior.alpha, prior.alphaKsi, prior.betaKsi}{Prior parameters
      for Gaussian mixture model that is calculated for each
      subnetwork (normal-inverse-Gamma prior). alpha tunes the mean;
      alphaKsi and betaKsi are the shape and scale parameters of the
      inverse Gamma function, respectively.}

  \item{vdp.threshold}{Minimal free energy improvement after which the
    variational Gaussian mixture algorithm is deemed converged.}

  \item{merging.threshold}{Minimal cost value improvement required for
  merging two subnetworks.}  

  \item{ite}{Defines maximum number of iterations on posterior
    update (updatePosterior). Increasing this can potentially lead to
    more accurate results, but computation may take longer.}

  \item{information.criterion}{Information criterion for model
  selection. Default is BIC (Bayesian Information Criterion); other
  options include AIC and AICc.}

  \item{speedup}{Takes advantage of approximations to PCA, mutual
    information etc in various places to speed up
    calculations. Particularly useful with large and densely connected
    networks and/or large sample size.}

  \item{speedup.max.edges}{Used if speedup = TRUE. Applies prefiltering
  of edges for calculating new joint models between subnetwork pairs
  when potential cost changes (delta) are updated for a newly merged
  subnetwork and its neighborghs. Empirical mutual information between
  each such subnetwork pair is calculated based on their first principal
  components, and joint models will be calculated only for the top
  candidates up to the number specified by speedup.max.edges. It is
  expected that the subnetwork pair that will benefit most from joint
  modeling will be among the top mutual infomation candidates. This way
  it is possible to avoid calculating exhaustive many models on the
  network hubs.}

  \item{mc.cores}{Number of cores to be used in parallelization. See
help(mclapply) for details.}
 
  \item{...}{Further optional arguments to be passed.}

}
  
\examples{

library(netresponse)
data( toydata )        # Load toy data set
D    <- toydata$emat   # Response matrix (for example, gene expression)
netw <- toydata$netw   # Network

# Run NetReponse algorithm
model <- detect.responses(D, netw, verbose = FALSE)

}
\keyword{ methods }
\keyword{ iteration }

\name{vdp.mixt}
\alias{vdp.mixt}
\title{vdp.mixt}
\description{Accelerated variational Dirichlet process Gaussian mixture.}
\usage{
vdp.mixt(dat, prior.alpha = 1, prior.alphaKsi = 0.01, prior.betaKsi = 0.01, do.sort = TRUE, threshold = 1e-05, initial.K = 1, ite = Inf, implicit.noise = 0, c.max = 10, speedup = TRUE)
}
\arguments{
    \item{dat}{Data matrix (samples x features).}

    \item{prior.alpha, prior.alphaKsi, prior.betaKsi}{Prior
      parameters for Gaussian mixture model (normal-inverse-Gamma
      prior). alpha tunes the mean; alphaKsi and betaKsi are the shape
      and scale parameters of the inverse Gamma function,
      respectively.}

    \item{do.sort}{When true, qOFz will be sorted in decreasing
      fashion by component size, based on colSums(qOFz). The qOFz
      matrix describes the sample-component assigments in the mixture
      model.}

    \item{threshold}{Defines the minimal free energy improvement that
    stops the algorithm: used to define convergence limit.}

    \item{initial.K}{Initial number of mixture components.}

    \item{ite}{Defines maximum number of iterations on posterior
    update (updatePosterior). Increasing this can potentially lead to
    more accurate results, but computation may take longer.}

    \item{implicit.noise}{Adds implicit noise; used by
    vdp.mk.log.lambda.so and vdp.mk.hp.posterior.so. By adding noise
    (positive values), one can avoid overfitting to local optima in
    some cases, if this happens to be a problem.}

    \item{c.max}{Maximum number of candidates to consider in
      find.best.splitting. During mixture model calculations new
      mixture components can be created until this upper limit has
      been reached. Defines the level of truncation for a truncated
      stick-breaking process.}

   \item{speedup}{When learning the number of components, each
	          component is splitted based on its first PCA
	          component. To speed up, approximate by using only
	          subset of data to calculate PCA.}

  }

\value{
  \item{ prior }{Prior parameters of the vdp-gm model.}
  \item{ posterior }{Posterior estimates for the model parameters and statistics. 

    weights: Mixture proportions, or weights, for the Gaussian mixture components.

    centroids: Centroids of the mixture components.

    sds: Standard deviations for the mixture model components (posterior modes of the covariance diagonals square root). Calculated as sqrt(invgam.scale/(invgam.shape + 1)).

    qOFz: Sample-to-cluster assigments (soft probabilistic associations).

    Nc: Component sizes    

    invgam.shape: Shape parameter (alpha) of the inverse Gamma distribution

    invgam.scale: Scale parameter (beta) of the inverse Gamma distribution

    Nparams: Number of model parameters

    K: Number of components in the mixture model
  }
  \item{ opts }{Model parameters that were used.}
  \item{ free.energy }{Free energy of the model.}
}

\details{Implementation of the Accelerated variational Dirichlet
  process Gaussian mixture model algorithm by Kenichi Kurihara et al.,
  2007.}

\references{
  Kenichi Kurihara, Max Welling and Nikos Vlassis: Accelerated
  Variational Dirichlet Process Mixtures.  In B. Sch\"olkopf and
  J. Platt and T. Hoffman (eds.), Advances in Neural Information
  Processing Systems 19, 761--768. MIT Press, Cambridge, MA 2007.
}
\author{Maintainer: Leo Lahti \email{leo.lahti@iki.fi}}

\note{This implementation is based on the Variational Dirichlet Process
  Gaussian Mixture Model implementation, Copyright (C) 2007 Kenichi
  Kurihara (all rights reserved) and the Agglomerative Independent
  Variable Group Analysis package (in Matlab): Copyright (C) 2001-2007
  Esa Alhoniemi, Antti Honkela, Krista Lagus, Jeremias Seppa, Harri
  Valpola, and Paul Wagner.}

\examples{
  set.seed(123)

  # Generate toy data with two Gaussian components
  dat <- rbind(array(rnorm(400), dim = c(200,2)) + 5,
               array(rnorm(400), dim = c(200,2)))

  # Infinite Gaussian mixture model with 
  # Variational Dirichlet Process approximation
  mixt <- vdp.mixt( dat )

  # Centroids of the detected Gaussian components
  mixt$posterior$centroids

  # Hard mixture component assignments for the samples
  apply(mixt$posterior$qOFz, 1, which.max)

}
\keyword{methods}
\keyword{iteration}

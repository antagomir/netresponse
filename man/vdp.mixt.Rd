\name{vdp.mixt}
\alias{vdp.mixt}
\title{vdp.mixt}
\description{Accelerated variational Dirichlet process Gaussian mixture.}
\usage{
vdp.mixt(dat, prior.alpha = 1, prior.alphaKsi = 0.01, prior.betaKsi = 0.01, do.sort = TRUE, threshold = 1e-05, initial.K = 1, ite = Inf, implicit.noise = 0, c.max = 10)
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
    vdp.mk.log.lambda.so and vdp.mk.hp.posterior.so By adding noise
    (positive values), one can avoid overfitting to local optima in
    some cases, if this happens to be a problem.}

    \item{c.max}{Maximum number of candidates to consider in
      find.best.splitting. During mixture model calculations new
      mixture components can be created until this upper limit has
      been reached. Defines the level of truncation for a truncated
      stick-breaking process.}
  }

\value{
  \item{ hp.prior }{Prior. }
  \item{ opts }{Optional model parameters that were used.}
  \item{ free.energy }{Free energy of the model.}
  \item{ hp.posterior }{Posterior estimates for the model parameters and statistics.}
  \item{ Kreal }{Number of non-empty mixture components.}
  \item{ qOFz }{Mixture component assignments for the samples;
  probability of each component, given sample.}
  \item{ Nparams }{Number of model parameters.}
  \item{ means }{Posterior modes for component centroids.}
  \item{ variances }{Posterior modes for component variances. Each
    component has a vector of variance parameters over the nodes, this
    determines the diagonal of the diagonal covariance matrices in the model)}
  \item{ weights }{Weight for each mixture component.}
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
\author{Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.
Maintainer: Leo Lahti <leo.lahti@iki.fi>}

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
  mixt$hp.posterior$Mubar

  # Mixture component assignments for the samples
  apply(mixt$hp.posterior$qOFz, 1, which.max)

}
\keyword{methods}
\keyword{iteration}

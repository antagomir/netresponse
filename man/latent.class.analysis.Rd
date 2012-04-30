\name{latent.class.analysis}
\alias{latent.class.analysis}
\title{Description: Latent class analysis based on (infinite) Gaussian mixture model. If the input (dat) is data matrix, a multivariate model is fitted. If the input is a vector or a 1-dimensional matrix, a univariate model is fitted.}
\usage{
  latent.class.analysis(dat, initial.responses = 1,
    max.responses = 10, verbose = TRUE, prior.alpha = 1,
    prior.alphaKsi = 0.01, prior.betaKsi = 0.01,
    update.hyperparams = 0, implicit.noise = 0,
    vdp.threshold = 1e-05, ite = Inf, speedup = TRUE,
    information.criterion = "BIC", ...)
}
\arguments{
  \item{dat}{vector (for univariate analysis) or a matrix
  (for multivariate analysis)}

  \item{initial.responses}{Initial number of components for
  each subnetwork model. Used to initialize calculations.}

  \item{max.responses}{Maximum number of responses for each
  subnetwork. Can be used to limit the potential number of
  network states.}

  \item{verbose}{Logical. Verbose parameter.}

  \item{implicit.noise}{Implicit noise parameter. Add
  implicit noise to vdp mixture model. Can help to avoid
  overfitting to local optima, if this appears to be a
  problem.}

  \item{update.hyperparams}{Logical. Indicate whether to
  update hyperparameters during modeling.}

  \item{prior.alpha}{Prior parameter for Gaussian mixture
  model. See help(detect.responses)}

  \item{prior.alphaKsi}{Prior parameter for Gaussian
  mixture model. See help(detect.responses)}

  \item{prior.betaKsi}{Prior parameter for Gaussian mixture
  model. See help(detect.responses)}

  \item{vdp.threshold}{Minimal free energy improvement
  after which the variational Gaussian mixture algorithm is
  deemed converged.}

  \item{ite}{Defines maximum number of iterations on
  posterior update (updatePosterior). Increasing this can
  potentially lead to more accurate results, but
  computation may take longer.}

  \item{information.criterion}{Information criterion for
  model selection. Default is BIC (Bayesian Information
  Criterion); other options include AIC and AICc.}

  \item{speedup}{Takes advantage of approximations to PCA,
  mutual information etc in various places to speed up
  calculations. Particularly useful with large and densely
  connected networks and/or large sample size.}

  \item{...}{Further optional arguments to be passed

  Returns:}
}
\value{
  Fitted latent class model (parameters and free energy)
}
\description{
  Arguments:
}
\author{
  Contact: Leo Lahti \email{leo.lahti@iki.fi}
}
\references{
  See citation("netresponse")
}
\keyword{utilities}


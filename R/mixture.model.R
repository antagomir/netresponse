# Save the other half of yourselves and your lives for pleasure and
# adventure. It is not enough to fight for natural land and the west;
# it is even more important to enjoy it. While you can. While it's
# still there... Enjoy yourselves, keep your brain in your head and
# your head firmly attached to the body, the body active and alive,
# and I promise you this much: I promise you this one sweet victory
# over our enemies, over those desk-bound men with their hearts in a
# safe deposit box, and their eyes hypnotized by desk calculators. I
# promise you this: you will outlive the bastards. --Ed Abbey.

#' @title Mixture model
#' @description Fit Gaussian mixture model
#' @param x data matrix (samples x features, for multivariate analysis) or a vector (for univariate analysis)
#' @param mixture.method Specify the approach to use in mixture modeling. Options. vdp (nonparametric Variational Dirichlet process mixture model); bic (based on Gaussian mixture modeling with EM, using BIC to select the optimal number of components)
#' @param max.responses Maximum number of responses for each subnetwork. Can be used to limit the potential number of network states.
#' @param implicit.noise Implicit noise parameter. Add implicit noise to vdp
#'        mixture model. Can help to avoid overfitting to local optima, if this
#'        appears to be a problem.
#' @param prior.alpha,prior.alphaKsi,prior.betaKsi Prior parameters for
#'        Gaussian mixture model that is calculated for each subnetwork
#'        (normal-inverse-Gamma prior). alpha tunes the mean; alphaKsi and betaKsi are
#'        the shape and scale parameters of the inverse Gamma function, respectively.
#' @param vdp.threshold Minimal free energy improvement after which the
#'        variational Gaussian mixture algorithm is deemed converged.
#' @param initial.responses Initial number of components for each subnetwork
#'         model. Used to initialize calculations.
#' @param ite Maximum number of iterations on posterior update
#'        (updatePosterior). Increasing this can potentially lead to more accurate results, but computation may take longer.
#' @param speedup Takes advantage of approximations to PCA, mutual information
#'         etc in various places to speed up calculations. Particularly useful with
#' 	   large and densely connected networks and/or large sample size.
#' @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture with mixture.method = "bic"
#' @param pca.basis pca.basis
#' @param min.responses minimum number of responses
#' @param ... Further optional arguments to be passed.
#' @return List with two elements: model: fitted mixture model (parameters and free energy); model.params: model parameters
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
mixture.model <- function (x, mixture.method = "vdp", max.responses = 10, implicit.noise = 0, prior.alpha = 1, prior.alphaKsi = 0.01, prior.betaKsi = 0.01, vdp.threshold = 1.0e-5, initial.responses = 1, ite = Inf, speedup = TRUE, bic.threshold = 0, pca.basis = FALSE, min.responses = 1, ...) {

  # Present data in PCA space to cope with diagonality of the covariances
  if (pca.basis) {    
    if (nrow(x) > ncol(x)) {
      x <- princomp(x)$scores
    } else {
      # FIXME: implement sparse PCA here to gain more generality?
      warning("Less samples than features, not applying PCA basis")
      x <- x
    }
  }

  if ( is.vector(x) ) { xm <- matrix(x, nrow = length(x)); rownames(xm) <- names(x); x <- xm }

  if (mixture.method == "vdp") {

    model <- vdp.mixt(x,
                      implicit.noise = implicit.noise,
		      prior.alpha = prior.alpha,
                      prior.alphaKsi = prior.alphaKsi,
		      prior.betaKsi = prior.betaKsi,
                      threshold = vdp.threshold,
		      initial.K = initial.responses, # FIXME: move initial.K into initial.responses for clarity
		      ite = ite,
		      c.max = max.responses,
		      speedup = speedup)

    model.params <- pick.model.parameters(model, colnames(x))      

  } else if (mixture.method == "bic") { 	

    model <- bic.mixture(x, max.modes = max.responses, bic.threshold = bic.threshold, min.modes = min.responses)  

    mu <- matrix(model$means, nrow = length(model$ws))
    sd <- matrix(model$sds, nrow = length(model$ws))
    ws <- matrix(model$ws)
    qofz <- matrix(model$qofz, ncol = length(ws))
    bic <- model$bic

    rownames(qofz) <- rownames(x)

    rownames(mu) <- rownames(sd) <- names(ws)
    colnames(mu) <- colnames(sd) <- colnames(x)        
   
    model.params <- list(mu = mu,
    		         sd = sd,
			 w = ws, 
			 qofz = qofz,
			 free.energy = model$free.energy, 
			 Nparams = model$Nparams, 
			 bic = bic)

  } else {
    stop("Provide proper mixture.method argument.")
  }    

  # FIXME: perhaps the 'model' is not needed any more when model.params is given? Check and remove to save space.
  # How about pca.basis = TRUE case?
  model.params

}


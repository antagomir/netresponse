# Copyright (C) 2008-2012 Olli-Pekka Huovilainen and Leo Lahti 
# Contact: Leo Lahti <leo.lahti@iki.fi>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

# Save the other half of yourselves and your lives for pleasure and
# adventure. It is not enough to fight for natural land and the west;
# it is even more important to enjoy it. While you can. While it's
# still there... Enjoy yourselves, keep your brain in your head and
# your head firmly attached to the body, the body active and alive,
# and I promise you this much: I promise you this one sweet victory
# over our enemies, over those desk-bound men with their hearts in a
# safe deposit box, and their eyes hypnotized by desk calculators. I
# promise you this: you will outlive the bastards. --Ed Abbey.

#' Description: Fit Gaussian mixture model
#'
#' Arguments:
#'  @param x data matrix (for multivariate analysis) or a vector (for univariate analysis)
#'  @param mixture.method Specify the approach to use in mixture modeling.
#'         Options. vdp (nonparametric Variational Dirichlet process mixture model);
#'         bic (based on Gaussian mixture modeling with EM, using BIC to select the
#'         optimal number of components)
#'  @param max.responses Maximum number of responses for each subnetwork. Can be
#'         used to limit the potential number of network states.
#'  @param implicit.noise Implicit noise parameter. Add implicit noise to vdp
#'        mixture model. Can help to avoid overfitting to local optima, if this
#'        appears to be a problem.
#'  @param prior.alpha,prior.alphaKsi,prior.betaKsi Prior parameters for
#'        Gaussian mixture model that is calculated for each subnetwork
#'        (normal-inverse-Gamma prior). alpha tunes the mean; alphaKsi and betaKsi are
#'        the shape and scale parameters of the inverse Gamma function, respectively.
#'  @param vdp.threshold Minimal free energy improvement after which the
#'        variational Gaussian mixture algorithm is deemed converged.
#'  @param initial.responses Initial number of components for each subnetwork
#'         model. Used to initialize calculations.
#'  @param ite Maximum number of iterations on posterior update
#'        (updatePosterior). Increasing this can potentially lead to more accurate results, but computation may take longer.
#'  @param speedup Takes advantage of approximations to PCA, mutual information
#'         etc in various places to speed up calculations. Particularly useful with
#' 	   large and densely connected networks and/or large sample size.
#'  @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture with mixture.method = "bic"
#'  @param ... Further optional arguments to be passed.
#'
#' Returns:
#'   @return List with two elements: model: fitted mixture model (parameters and free energy); model.params: model parameters
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

mixture.model <- function (x, mixture.method = "vdp", max.responses = 10, implicit.noise = 0, prior.alpha = 1, prior.alphaKsi = 0.01, prior.betaKsi = 0.01, vdp.threshold = 1.0e-5, initial.responses = 1, ite = Inf, speedup = TRUE, bic.threshold = 0, ...) {

  if (mixture.method == "vdp") {

    if ( is.vector(x) ) { x <- matrix(x, nrow = length(x)) }

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
    model <- bic.mixture(x, max.modes = max.responses, bic.threshold = bic.threshold)  
    model.params <- list(mu = model$means, sd = model$sds, w = model$ws, free.energy = model$free.energy, Nparams = model$Nparams)
  } else {
    stop("Provide proper mixture.method argument.")
  }    

  list(model = model, params = model.params)

}







#' Description: Latent class analysis based on (infinite) Gaussian mixture model. 
#' If the input is data matrix, a multivariate model is fitted; if the input is a vector, a univariate model is fitted
#'
#' Arguments:
#'  @param x samples x features matrix for multivariate analysis, or a vector for univariate analysis 
#'  @param max.modes Maximum number of modes to be checked for mixture model selection
#'  @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture.
#'  @param ... Further optional arguments to be passed
#'
#' Returns:
#'   @return Fitted latent class model (parameters and free energy)
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
bic.mixture <- function (x, max.modes, bic.threshold = 0, ...) { 

  if (is.vector(x)) {
    bic.mixture.univariate(x, max.modes, bic.threshold, ...)
  } else {
    bic.mixture.multivariate(x, max.modes, bic.threshold, ...)
  }

}


#' Description: Latent class analysis based on (infinite) Gaussian mixture model. 
#' If the input (dat) is data matrix, a multivariate model is fitted. 
#'
#' Arguments:
#'  @param x  matrix (for multivariate analysis)
#'  @param max.modes Maximum number of modes to be checked for mixture model selection
#'  @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture.
#'  @param ... Further optional arguments to be passed
#'
#' Returns:
#'   @return Fitted latent class model (parameters and free energy)
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

bic.mixture.multivariate <- function (x, max.modes, bic.threshold = 0, ...) { 

  #x <- mat; max.modes = params$max.responses; bic.threshold = params$bic.threshold

  best.mode <- bic.select.best.mode(x, max.modes, bic.threshold) 

  mcl <- Mclust(x, G = best.mode)

  means <- t(mcl$parameters$mean)
  vars <- t(apply(mcl$parameters$variance$sigma,3, function(x){diag(x)}))
  sds <- sqrt(vars)
  ws <- as.vector(mcl$parameters$pro)
  if (is.null(ws)) {ws <- 1} 

  Nparams <- prod(dim(means)) + prod(dim(sds)) + length(ws) 

  rownames(means) <- rownames(sds) <- names(ws) <- paste("Mode", 1:length(ws), sep = "-")
  colnames(means) <- colnames(sds) <- colnames(x)

  list(means = means, sds = sds, ws = ws, Nparams = Nparams, free.energy = -mcl$loglik)

}



#' Description: Select optimal number of mixture components by adding components until 
#' the increase in objective function is below threshold.
#'
#' Arguments:
#'  @param x  dat vector (for univariate analysis) or a matrix (for multivariate analysis)
#'  @param max.modes Maximum number of modes to be checked for mixture model selection
#'  @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture.
#'
#' Returns:
#' @return Fitted latent class model (parameters and free energy)
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
bic.select.best.mode <- function (x, max.modes, bic.threshold) {

  # Cost for single mode
  # BIC : smaller is better
  # mclustBIC returns the value for -BIC, to be exact
  nc <- 1
  if (is.vector(x)) { # univariate
    m <- -mclustBIC(x, G = nc)[, "V"] 
  } else { # multivariate
    m <- -mclustBIC(x, G = nc)[, "VVV"] # BIC : smaller is better
  }
  
  # ----------------------------------------------------------------
  
  add.component <- TRUE
  best.mode <- 1
  if (max.modes == 1) {
    add.component <- FALSE
  }

  while (add.component && nc < max.modes) {

    nc <- nc + 1

    if (is.vector(x)) { # univariate
      m.new <- -mclustBIC(x, G = nc)[, "V"] 
    } else { # multivariate
      m.new <- -mclustBIC(x, G = nc)[, "VVV"] # BIC : smaller is better
    }

    # FIXME: compressing data with PCA after dimensionality gets otherwise too high?
    # with around ncol(x) = 30 the mclustBIC is starting to produce NAs
    if (is.na(m.new)) {save(x, nc, file = "m.new.RData")}
    
    bic.delta <- m.new - m

    if (bic.delta < -bic.threshold) { 
      best.mode <- nc 
      m <- m.new
    } else {
      add.component <- FALSE
    }
  }

  best.mode

}



#' Description: Latent class analysis based on (infinite) Gaussian mixture
#' model. If the input (dat) is data matrix, a multivariate model is fitted. If
#' the input is a vector or a 1-dimensional matrix, a univariate model is
#' fitted.
#'
#' Arguments:
#'  @param x  dat vector (for univariate analysis) or a matrix (for multivariate analysis)
#'  @param max.modes Maximum number of modes to be checked for mixture model selection
#'  @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture.
#'  @param ... Further optional arguments to be passed
#'
#' Returns:
#' @return Fitted latent class model (parameters and free energy)
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
bic.mixture.univariate <- function (x, max.modes, bic.threshold = 0, ...) { 

  # x <- datamatrix[, node];  max.modes = params$max.responses; bic.threshold = params$bic.threshold

  library(mclust)

  best.mode <- bic.select.best.mode(x, max.modes, bic.threshold) 
  mcl <- Mclust(x, G = best.mode)

  means <- as.vector(mcl$parameters$mean)
  sds <- as.vector(sqrt(mcl$parameters$variance$sigmasq))
  if (length(sds) == 1) {sds <- rep(sds, length(means))} 
  ws <- as.vector(mcl$parameters$pro)
  if (is.null(ws)) {warning("NULL weights, replacing with 1"); ws <- 1} 
  if (is.null(means)) {warning("NULL means, replacing with 1"); means <- 1} 
  if (is.null(sds)) {warning("NULL sds, replacing with 1"); sds <- 1} 

  Nparams <- length(means) + length(sds) + length(ws) 

  names(means) <- names(sds) <- names(ws) <- paste("Mode", 1:length(ws), sep = "-")

  list(means = means, sds = sds, ws = ws, Nparams = Nparams, free.energy = -mcl$loglik)

}


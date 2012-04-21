# Copyright (C) 2012 Leo Lahti
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

#' Description: Latent class analysis based on (infinite) Gaussian mixture model. If the input (dat) is data matrix, a multivariate model is fitted. If the input is a vector or a 1-dimensional matrix, a univariate model is fitted.
#'
#' Arguments:
#'   @param dat vector (for univariate analysis) or a matrix (for multivariate analysis)
#'   @param initial.responses Initial number of components for each subnetwork model. Used to initialize calculations.
#'   @param max.responses  Maximum number of responses for each
#'          subnetwork. Can be used to limit the potential number of network
#'          states.
#'  @param verbose Logical. Verbose parameter.
#'  @param implicit.noise Implicit noise parameter. Add implicit noise
#'         to vdp mixture model. Can help to avoid overfitting to local optima,
#'         if this appears to be a problem.
#'  @param update.hyperparams Logical. Indicate whether to update hyperparameters during modeling.
#'  @param prior.alpha Prior parameter for Gaussian mixture model. See help(detect.responses)
#'  @param prior.alphaKsi Prior parameter for Gaussian mixture model. See help(detect.responses)
#'  @param prior.betaKsi Prior parameter for Gaussian mixture model. See help(detect.responses)
#'  @param vdp.threshold Minimal free energy improvement after which the variational Gaussian mixture algorithm is deemed converged.
#'  @param ite Defines maximum number of iterations on posterior
#'    	   update (updatePosterior). Increasing this can potentially lead to
#'    	   more accurate results, but computation may take longer.
#'  @param information.criterion Information criterion for model
#'  	   selection. Default is BIC (Bayesian Information Criterion); other
#'  	   options include AIC and AICc.
#'  @param speedup Takes advantage of approximations to PCA, mutual
#'    	   information etc in various places to speed up
#'    	   calculations. Particularly useful with large and densely connected
#'    	   networks and/or large sample size.
#'  @param ... Further optional arguments to be passed
#'
#'
#' Returns:
#'   @return Fitted latent class model (parameters and free energy)
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities

         
latent.class.analysis <- function(dat,
         initial.responses = 1,   # initial number of components. FIXME: is this used?
         max.responses = 10,      
         verbose = TRUE,          # print proc. information
         prior.alpha    = 1,      # Prior parameters
         prior.alphaKsi = 0.01,   # for VDP mixture
         prior.betaKsi  = 0.01,   # scale parameter for inverse Gamma
         update.hyperparams = 0,  # update hyperparameters. FIXME: check if this is applicable.
         implicit.noise = 0,      # Add implicit noise in vdp.mk.log.lambda.so and vdp.mk.hp.posterior.so 
         vdp.threshold = 1.0e-5,  # min. free energy improvement that stops VDP
         ite = Inf,               # max. iterations in updatePosterior
	 speedup = TRUE,
	 information.criterion = "BIC",
         ... # Further arguments
)

{

  # For debugging: dat <- rnorm(100); initial.responses = 1; max.responses = 10; max.subnet.size = 10; verbose = TRUE; prior.alpha = 1; prior.alphaKsi = 0.01; prior.betaKsi  = 0.01; update.hyperparams = 0; implicit.noise = 0; vdp.threshold = 1.0e-5; merging.threshold = 0; ite = Inf; information.criterion = "AIC"; speedup = TRUE

  # If input data is a vector, convert it into a matrix
  if (is.vector(dat)) { dat <- as.matrix(dat) } 

  # Check that the data matrix is of correct format
  datamatrix <- check.matrix(dat)  

  # Store all params 
  params <- list(initial.responses = initial.responses, 
  	    	 max.responses = max.responses,
                 verbose = verbose, 
		 prior.alpha = prior.alpha, 
		 prior.alphaKsi = prior.alphaKsi,
                 prior.betaKsi = prior.betaKsi, 
		 update.hyperparams = update.hyperparams, 
		 implicit.noise = implicit.noise,
                 vdp.threshold = vdp.threshold,
		 ite = ite, 
		 Nlog = log( nrow( datamatrix ) ),
		 speedup = speedup,
		 information.criterion = information.criterion
		 )

  # Fit Gaussian mixture model
  model <- vdp.mixt( datamatrix,
                      implicit.noise = params$implicit.noise,
                      prior.alpha = params$prior.alpha,
                      prior.alphaKsi = params$prior.alphaKsi,
                      prior.betaKsi = params$prior.betaKsi,
                      threshold = params$vdp.threshold,
                      initial.K = params$initial.responses,
                      ite = params$ite,
                      c.max = params$max.responses - 1,
                      speedup = params$speedup )

  # Cost for model, including penalty for model parameters for comparing models in a later netresponse merge stage
  # C <- info.criterion(model$posterior$Nparams, params$Nlog, -model$free.energy, criterion = params$information.criterion) 

  # Model 
  fitted.model <- pick.model.parameters(model, nodes = colnames(datamatrix))

  fitted.model

}

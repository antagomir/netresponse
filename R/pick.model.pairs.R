# Copyright (C) 2008-2012 Leo Lahti and Olli-Pekka Huovilainen
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
# Acknowledgements: This program is based on the AIVGA Agglomerative
# Independent Variable Group Analysis package (v. 1.0) Copyright (C)
# 2001-2007 Esa Alhoniemi, Antti Honkela, Krista Lagus, Jeremias
# Seppa, Harri Valpola, and Paul Wagner.


#' @title Pick model pairs
#' @description Mainly for internal use. Calculate joint model for each node pair
#' @param network network
#' @param network.nodes network.nodes
#' @param node.models node.models
#' @param datamatrix datamatrix
#' @param params parameters
#' @return \item{model.pairs}{joint models for each node pair}
#' \item{delta}{corresponding delta value for the cost function}
#' @author Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @export
#' @examples #
pick.model.pairs <- function (network, network.nodes, node.models, datamatrix, params) {

  # Storage list for calculated models
  model.pairs <- vector(length = ncol(network), mode = "list" ) 
  delta <- rep(NA, ncol(network))

  for (edge in 1:ncol(network)){

    tmp <- edge.delta(edge, network = network, 
    	   		    network.nodes = network.nodes, 
	   		    datamatrix = datamatrix, 
			    params = params,
	   		    node.models = node.models)

    model.pairs[[edge]] <- tmp[["model"]]
    delta[[edge]] <- tmp[["delt"]]

  }

  gc()

  list(model.pairs = model.pairs, delta = delta)

}



edge.delta <- function (edge, network, network.nodes, datamatrix, params, node.models) {

    if ( params$verbose ) { message(paste('Computing delta values for edge ', edge, '/', ncol(network), '\n')) }
    a <- network[1, edge]
    b <- network[2, edge]
    vars <- network.nodes[c(a, b)]

    mat <- matrix(datamatrix[, vars], nrow( datamatrix ))
    rownames(mat) <- rownames(datamatrix)
    colnames(mat) <- vars

    model.params <- mixture.model(x = mat, 
    	   		 mixture.method = params$mixture.method, 
			 max.responses = params$max.responses, 
			 implicit.noise = params$implicit.noise, 
			 prior.alpha = params$prior.alpha, 
			 prior.alphaKsi = params$prior.alphaKsi, 
			 prior.betaKsi = params$prior.betaKsi, 
			 vdp.threshold = params$vdp.threshold, 
			 initial.responses = params$initial.responses, 
			 ite = params$ite, 
			 speeup = params$speedup,
			 bic.threshold = params$bic.threshold, 
			 pca.basis = params$pca.basis
			 ) 

    # Compute COST-value for two independent subnets vs. joint model
    # Negative free energy (-cost) is (variational) lower bound for P(D|H)
    # Use it as an approximation for P(D|H)
    # Cost for the indpendent and joint models
    # -cost is sum of two independent models (cost: appr. log-likelihoods)
    costind     <-  info.criterion(node.models[[a]]$Nparams + node.models[[b]]$Nparams, params$Nlog, -(node.models[[a]]$free.energy + node.models[[b]]$free.energy), criterion = params$information.criterion)
    costjoint   <-  info.criterion(model.params$Nparams, params$Nlog, -model.params$free.energy, criterion = params$information.criterion)

    # NOTE: COST is additive so summing is ok
    # change (increase) of the total COST / cost
    delt <- as.numeric(costjoint - costind)

    # Store these only if it would improve the cost; otherwise never needed again
    if (-delt > params$merging.threshold) {
      mod.pair <- model.params
    } else {
      mod.pair <- 0
    }
			
    return(list(model = mod.pair, delt = delt))
}




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
# Acknowledgements: This program is based on the AIVGA Agglomerative
# Independent Variable Group Analysis package (v. 1.0) Copyright (C)
# 2001-2007 Esa Alhoniemi, Antti Honkela, Krista Lagus, Jeremias
# Seppa, Harri Valpola, and Paul Wagner.



#' independent.models
#' 
#' Mainly for internal use. Provide independent models for each node.
#'
#' @param datamatrix datamatrix
#' @param params parameters
#' @return 
#'   \item{nodes }{Model for each node} 
#'   \item{C }{Costs for individual models}
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @examples #
#' @keywords internal
independent.models <- function (datamatrix, params) {

  # Storage list for calculated models
  model.nodes <- vector(length = ncol(datamatrix), mode = "list" ) # individual nodes

  if (params$verbose) { message("Compute cost for each variable") }

  C <- vector(length = ncol(datamatrix), mode = "numeric")

  # FIXME parallelize?
  for (k in 1:ncol(datamatrix)){

    node <- colnames(datamatrix)[[k]]
    
    if ( params$verbose ) { message(paste('Computing model for node', k, "/", ncol( datamatrix ))) }
    
    Nparams <- NULL
    if (params$mixture.method == "vdp") {

      model <- vdp.mixt( matrix(datamatrix[, node], nrow( datamatrix )),
                      implicit.noise = params$implicit.noise,
                      prior.alpha = params$prior.alpha,
                      prior.alphaKsi = params$prior.alphaKsi,
                      prior.betaKsi = params$prior.betaKsi,
                      threshold = params$vdp.threshold,
                      initial.K = params$initial.responses,
                      ite = params$ite,
                      c.max = params$max.responses - 1,
                      speedup = params$speedup )

      model.params <- pick.model.parameters(model, node)

    } else if (params$mixture.method == "bic") {

      model <- bic.mixture.univariate(datamatrix[, node], max.modes = params$max.responses, bic.threshold = params$bic.threshold)

      means <- matrix(model$means, nrow = length(model$ws))   
      sds <- matrix(model$sds, nrow = length(model$ws))   
      ws <- model$ws

      rownames(means) <- rownames(sds) <- names(ws) <- paste("Mode", 1:length(ws), sep = "-")
      colnames(means) <- colnames(sds) <- node
      
      model.params <- list(mu = means, sd = sds, w = ws, free.energy = model$free.energy, Nparams = model$Nparams)

    } else {
      stop("Provide proper mixture.method")
    }

    # Cost for model
    C[[k]] <- info.criterion(model.params$Nparams, 
    	      		params$Nlog, 
			-model.params$free.energy, 
			criterion = params$information.criterion) 

    model.nodes[[k]] <- model.params

  }

  gc()

  if (params$verbose) { message('independent models done') }

  list(nodes = model.nodes, C = C)  

}




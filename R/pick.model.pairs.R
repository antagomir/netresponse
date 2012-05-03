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


#' pick.model.pairs
#' 
#' Mainly for internal use. Calculate joint model for each node pair
#' 
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
	   		    datamatrix = datamatrix, params = params,
	   		    node.models = node.models)

    model.pairs[[edge]] <- tmp[["model"]]
    delta[[edge]] <- tmp[["delt"]]

  }

  gc()

  list(model.pairs = model.pairs, delta = delta)

}



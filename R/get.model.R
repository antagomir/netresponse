#' get.model
#' 
#' Mainly for internal use. Calculate joint model for a given node pair.
#' 
#' 
#' @usage get.model(datamatrix, network, edge, network.nodes, G, params)
#' @param datamatrix datamatrix
#' @param network network
#' @param edge edge
#' @param network.nodes network.nodes
#' @param G G
#' @param params parameters
#' @return NetResponse model object
#' @author Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @examples #
get.model <- function (datamatrix, network, edge, network.nodes, G, params) { 
	
  # Pick node indices          
  a <- network[1, edge]          
  i <- network[2, edge]           
  vars  <- network.nodes[sort(c(G[[a]], G[[i]]))]          

  # model <- 
  vdp.mixt(matrix(datamatrix[, vars], nrow( datamatrix )),          
                          implicit.noise = 0,
                          prior.alpha = params$prior.alpha,
                          prior.alphaKsi = params$prior.alphaKsi,
                          prior.betaKsi = params$prior.betaKsi,
                          threshold = params$vdp.threshold,
                          initial.K = params$initial.responses,
                          ite = params$ite,
                          c.max = params$max.responses,
                          speedup = params$speedup)
  
}
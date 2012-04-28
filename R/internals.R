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




#' check.network
#' 
#' Mainly for internal use. Check input network and make formatting for
#' detect.responses
#' 
#' 
#' @usage check.network(network, datamatrix, verbose = FALSE)
#' @param network Input network, see detect.responses
#' @param datamatrix Input datamatrix, see detect.responses
#' @param verbose Print intermediate messages
#' @return \item{formatted }{Formatted network (self-links removed)}
#' \item{original }{Original network (possible in another representation
#' format)} \item{delta }{Cost function changes corresponding to the
#' 'formatted' network.} \item{nodes }{Nodes corresponding to the 'formatted'
#' network.}
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.  Maintainer:
#' Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso detect.responses
#' @references See citation("netresponse")
#' @keywords internal
#' @examples
#' 
#' # TBA
#' 
check.network <- function (network, datamatrix, verbose = FALSE) {

  # If no network is given, assume fully connected net
  if (is.null(network)) { 
    if (verbose) { warning("No network provided in function call: assuming fully connected nodes.") }
    network <- matrix(1, ncol(datamatrix), ncol(datamatrix)) 
    rownames(network) <- colnames(network) <- colnames(datamatrix)
  }

  # FIXME: later add other forms of sparse matrices from Matrix package   
  accepted.formats.net <- c("matrix", "Matrix", "dgCMatrix", "dgeMatrix", "graphNEL", "igraph", "graphAM")
  if (!class(network)[[1]] %in% accepted.formats.net) {  
    stop(paste("network needs to be in one of the following formats:", paste(accepted.formats.net, collapse = "; ")))
  }
  
  # Convert matrix into graphNEL if it is not already
  if (!class(network) == "graphNEL") {
    if (class(network) %in% c("dgeMatrix", "dfCMatrix", "Matrix", "data.frame")) {
      # Sparse matrix or data.frame needs to be converted first into matrix 
      network <- as.matrix(network)
    }
    
    if (is.matrix(network)) {

      if ( !nrow(network) == ncol(network) ) { stop("Error: network nrow = ncol required.\n") }

      # Ensuring symmetric network
      if (any(!network == t(network))) {
        warning("Network is not symmetric. Removing link directions to force symmetric network.")
        network <- ((network + t(network)) > 0) - 0
      }
      
      # check that node names given in the data and correspond
      if ( is.null(rownames( network )) || is.null(colnames( datamatrix )) ) {        
        if ( !nrow(network) == ncol(datamatrix) ) {
          stop("Error: Equal number of features required for the network and data matrix when feature names not given.\n")
        } else {
          warning("Warning: network and/or data features are not named; matched by order.\n")
          if (is.null(rownames( network )) && is.null(colnames( datamatrix ))) {
            rownames(network) <- colnames(network) <- as.character(1:nrow(network))
            colnames(datamatrix) <- rownames(network)
          } else if (is.null(rownames( network )) && !is.null(colnames( datamatrix ))) {
            rownames(network) <- colnames(network) <- colnames( datamatrix )
          } else if (!is.null(rownames( network )) && is.null(colnames( datamatrix ))) {
            colnames( datamatrix ) <- rownames(network)
          }
        }
      }
      network <- as(new("graphAM", adjMat = network), "graphNEL")
    } else if (class(network) == "igraph") {
      network <- igraph.to.graphNEL(network)
    }
  }

  # Now the network is in graphNEL format

  # FIXME: adjust such that igraph does not need to be converted in graphNEL (which is larger
  # FIXME: add option to give this as input; seems to consume much less memory than graphNEL  

  # list network nodes that are not in datamatrix  
  common.feats <- intersect(nodes(network), colnames(datamatrix))
  other.feats <- setdiff(nodes(network), common.feats)

  # remove features with no functional data
  if (length(other.feats) > 0) {
    if (verbose) { message(paste("removing network nodes that are not in datamatrix: ", length(other.feats), " nodes removed;", length(common.feats), " nodes used for modeling.")) }
    #message("converting to igraph")
    network <- igraph.from.graphNEL(network)
    #message("selecting nodes that have functional data")
    network <- subgraph(network, common.feats)
    #message("converting to graphNEL")
    network <- igraph.to.graphNEL(network)
    #network <- removeNode(other.feats, network)    
  }

  message("convert the network into edge matrix")
  # store original network node list
  network.nodes <- nodes(network)
  network <- edgeMatrix(network, duplicates = FALSE) # indices correspond to node list in network.nodes
  # order such that row1 < row2
  network <- apply(network, 2, sort)
  if (verbose) message("removing self-links")  
  network <- network[, !network[1,] == network[2,]]

  # Store the network in igraph format  
  network.orig <- network # store network used for modeling (preprocessed)
  # FIXME: igraph is more memory-efficient but could not be used as network class in NetResponseModel definition
  # for some reason. If possible, convert from graphNEL to igraph later on.

  # Network rows:                                        
  # 1) nodes 2) nodes 3) merging cost function delta
  rownames(network) <- c("node1", "node2")

  # Delta corresponds to the columns of the 'network' object
  delta <- rep(NA, ncol(network))

  list(formatted = network, original = network.orig, delta = delta, nodes = network.nodes)
}





#' independent.models
#' 
#' Mainly for internal use. Provide independent models for each node.
#' 
#' 
#' @usage independent.models(datamatrix, params, mixture.method = "vdp")
#' @param datamatrix datamatrix
#' @param params parameters
#' @param mixture.method Specify the approach to use in mixture modeling.
#' Options. vdp (nonparametric Variational Dirichlet process mixture model);
#' bic (based on Gaussian mixture modeling with EM, using BIC to select the
#' optimal number of components)
#' @return \item{nodes }{Model for each node} \item{C }{Costs for individual
#' models}
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @export
#' @examples
independent.models <- function (datamatrix, params, mixture.method = "vdp") {

  # Storage list for calculated models
  model.nodes <- vector(length = ncol(datamatrix), mode = "list" ) # individual nodes
  
  if (params$verbose) { message("Compute cost for each variable") }

  C <- vector(length = ncol(datamatrix), mode = "numeric")

  # FIXME parallelize?
  for (k in 1:ncol(datamatrix)){

    node <- colnames(datamatrix)[[k]]
    
    if ( params$verbose ) { message(paste('Computing model for node', k, "/", ncol( datamatrix ))) }
    
    Nparams <- NULL
    if (mixture.method == "vdp") {

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

    } else if (mixture.method == "bic") {
      # FIXME: add c.max here

      model <- bic.mixture.univariate(datamatrix[, node], max.modes = params$max.responses)

      means <- matrix(model$means, nrow = length(model$means))   
      sds <- matrix(model$sds, nrow = length(model$sds))   
      ws <- matrix(model$ws, nrow = length(model$ws))   

      rownames(means) <- rownames(sds) <- names(ws) <- paste("Mode", 1:length(ws), sep = "-")
      colnames(means) <- colnames(sds) <- node
      
      model.params <- list(mu = means, sd = sds, w = ws, free.energy = model$free.energy, Nparams = model$Nparams)

    } else {
      stop("Provide proper mixture.method")
    }

    # COST for model

    C[[k]] <- info.criterion(model.params$Nparams, params$Nlog, -model.params$free.energy, criterion = params$information.criterion) 

    model.nodes[[k]] <- model.params

  }
  
  gc()

  if (params$verbose) { message('done') }

  list(nodes = model.nodes, C = C)  

}



#' pick.model.pairs
#' 
#' Mainly for internal use. Calculate joint models for each node pair
#' 
#' 
#' @usage pick.model.pairs(network, network.nodes, model.nodes, datamatrix, params, mixture.method = "vdp")
#' @param network network
#' @param network.nodes network.nodes
#' @param model.nodes model.nodes
#' @param datamatrix datamatrix
#' @param params parameters
#' @param mixture.method Specify the approach to use in mixture modeling.
#' Options. vdp (nonparametric Variational Dirichlet process mixture model);
#' bic (based on Gaussian mixture modeling with EM, using BIC to select the
#' optimal number of components)
#' @return \item{model.pairs}{joint models for each node pair}
#' \item{delta}{corresponding delta value for the cost function}
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.  Maintainer:
#' Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @export
#' @examples
pick.model.pairs <- function (network, network.nodes, model.nodes, datamatrix, params, mixture.method = "vdp") {

  # Storage list for calculated models
  model.pairs <- vector(length = ncol(network), mode = "list" ) 
  delta <- rep(NA, ncol(network))

  if (params$max.subnet.size > 1) {

    if (params$mc.cores == 1) {

      for (edge in 1:ncol(network)){

        if (params$verbose) { 
	  message(paste('Computing delta values for edge ', edge, '/', ncol(network), '\n')) 
	}

        a <- network[1, edge]
        b <- network[2, edge]  
        vars  <- network.nodes[c(a, b)]

        tmp <- mixture.model(matrix(datamatrix[, vars], nrow( datamatrix )), vars, params)
  	model <- tmp$model # FIXME: perhaps the 'model' is not needed when model.params is given. Check and remove.
  	model.params <- tmp$params

        # Compute COST-value for two independent subnets vs. joint model 
        # Negative free energy (-cost) is (variational) lower bound for P(D|H)
        # Use it as an approximation for P(D|H)
        # Cost for the indpendent and joint models
        # -cost is sum of two independent models (cost: appr. log-likelihoods)
        costind.ab <- info.criterion(model.nodes[[a]]$Nparams + model.nodes[[b]]$Nparams, params$Nlog, -(model.nodes[[a]]$free.energy + model.nodes[[b]]$free.energy), criterion = params$information.criterion)
        costjoint.ab <- info.criterion(model.params$Nparams, params$Nlog, -model.params$free.energy, criterion = params$information.criterion)

        # NOTE: COST is additive -> summing the scores is ok
        # change (increase) of the total COST / cost
        delta <- as.numeric(costjoint.ab - costind.ab)
      
        # Store these only if it would improve the cost; otherwise never needed
        if (-delta > params$merging.threshold) {    
          model.pairs[[edge]] <- model.params
        } else {
          model.pairs[[edge]] <- 0
        }
      }
    
    } else {

      if (params$verbose) {
        message(paste('Computing delta values for edges with multiple cores\n'))
        message(paste("Using", params$mc.cores, "cores", "\n")) 
      }
   
      # FIXME: parallelize - mclapply caused some problems with test/test.R		  
      res <- lapply(1:ncol(network), function (edge) {	
    	edge.delta(edge, network = network, network.nodes = network.nodes, 
			 datamatrix = datamatrix, params = params,
			 model.nodes = model.nodes, model = model)[[1]]})

      # Convert to vector
      delta <- unlist(lapply(res, function (x) {x$delt})) 

      # Convert to list
      model.pairs <- lapply(res, function (x) {x[1:5]}) # FIXME: parallelize

    }
  }

  gc()

  list(model.pairs = model.pairs, delta = delta)

}


#' update.model.pair
#' 
#' Mainly for internal use. Calculate joint model for given node pair and
#' update delta accordingly.
#' 
#' 
#' @usage update.model.pair(datamatrix, delta, network, edge, network.nodes, G,
#' params, model.nodes, model.pairs)
#' @param datamatrix datamatrix
#' @param delta delta
#' @param network network
#' @param edge edge
#' @param network.nodes network.nodes
#' @param G G
#' @param params params
#' @param model.nodes model.nodes
#' @param model.pairs model.pairs
#' @return \item{model.pairs }{model.pairs} \item{delta }{delta}
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.  Maintainer:
#' Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @examples
update.model.pair <- function (datamatrix, delta, network, edge, network.nodes, G, params, model.nodes, model.pairs) {

  # Pick node indices          
  a <- network[1, edge]          
  b <- network[2, edge]           
  vars  <- network.nodes[sort(c(G[[a]], G[[b]]))]          

  tmp <- mixture.model(matrix(datamatrix[, vars], nrow( datamatrix )), vars, params) 
  model <- tmp$model # FIXME: perhaps the 'model' is not needed when model.params is given. Check and remove.
  model.params <- tmp$params

  # Negative free energy is (variational) lower bound for P(D|H)          
  # Use this to approximate P(D|H)          
  if (is.finite(model$free.energy)) {
    # Compute COST-value for two independent subnets vs. joint model
    # Negative free energy (-cost) is (variational) lower bound for P(D|H)
    # Use it as an approximation for P(D|H)
    # Cost for the indpendent and joint models
    # -cost is sum of two independent models (cost: appr. log-likelihoods)
    cost.ind     <-  info.criterion(model.nodes[[a]]$Nparams + model.nodes[[b]]$Nparams, params$Nlog, -(model.nodes[[a]]$free.energy + model.nodes[[b]]$free.energy), criterion = params$information.criterion)
    cost.joint   <-  info.criterion(model.params$Nparams, params$Nlog, -model.params$free.energy, criterion = params$information.criterion)
    # change (increase) of the total cost
    delta[[edge]] <- cost.joint - cost.ind             
  } else  {
    warning("No free energy obtained.")            
    delta[[edge]] <- Inf       
  }          
          
  # Store the joint models / cost for two independent vs. joint model  
  if (-delta[[edge]] > params$merging.threshold) {  
    # Store joint model only if it would improve the cost            
    model.pairs[[edge]] <- model.params
  } else {          
    model.pairs[[edge]] <- 0
  }

  list(model.pairs = model.pairs, delta = delta)
  
}





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
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.  Maintainer:
#' Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @examples
#' 
#' # TBA
#' 
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
                          c.max = params$max.responses - 1,
                          speedup = params$speedup)
  
}



#' get.mis
#' 
#' Mainly for internal use. Estimate mutual information for node pairs based on
#' the first principal components
#' 
#' 
#' @usage get.mis(datamatrix, network, delta, network.nodes, G, params)
#' @param datamatrix datamatrix
#' @param network network
#' @param delta delta
#' @param network.nodes network.nodes
#' @param G G
#' @param params params
#' @return mutual information matrix
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.  Maintainer:
#' Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @examples
#' 
#' # TBA
#' 
get.mis <- function (datamatrix, network, delta, network.nodes, G, params) {

  require(minet)          
  mis <- c()          
  mi.cnt <- 0            
  for (edge in which(is.na(delta))){          
    mi.cnt <- mi.cnt + 1             
    # Pick node indices            
    a <- network[1, edge]            
    i <- network[2, edge]            
    dat <- cbind(prcomp(matrix(datamatrix[, network.nodes[G[[a]]]], nrow(datamatrix)), center = TRUE)$x[, 1],
                 prcomp(matrix(datamatrix[, network.nodes[G[[i]]]], nrow(datamatrix)), center = TRUE)$x[, 1])

    mis[[mi.cnt]] <- build.mim(dat, estimator="mi.empirical", disc = "equalwidth", nbins = params$nbins)[1, 2]

  }

  mis
}





#' filter.netw
#' 
#' Mostly for internal use. Prefilter edges if speedups required.
#' 
#' Include to the network only the edges with the highest mutual information,
#' calculated based on the first principal components.
#' 
#' @usage filter.netw(network, delta, datamatrix, params)
#' @param network network
#' @param delta associated cost function value changes for each node merge
#' @param datamatrix datamatrix
#' @param params parameters
#' @return Filtered network
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.  Maintainer:
#' Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @examples
#' 
#' # TBA
#' 
filter.netw <- function (network, delta, datamatrix, params) {

  tmp <- list(network = network, delta = delta)

  if (params$max.subnet.size > 1) {
    if (params$verbose) {message("Filter the network to only keep the edges with highest mutual information")}
    # Filter out the least promising edges from the network
    # based on mutual information. For each variable, pick at most speedup.max.edges
    if (params$speedup && !is.null(params$speedup.max.edges)) {
      tmp <- filter.network(network, delta, datamatrix, params)
    } 
  }

  tmp

}





#' check.matrix
#' 
#' Mostl for internal purposes. Check input matrix format.
#' 
#' 
#' @usage check.matrix(datamatrix)
#' @param datamatrix See detect.responses
#' @return The datamatrix, possibly added with necessary formatting for the
#' netresponse algorithm.
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.  Maintainer:
#' Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso detect.responses
#' @references See citation("netresponse")
#' @keywords internal
#' @examples
#' 
#' # datamatrix <- check.matrix(datamatrix)
#' 
check.matrix <- function (datamatrix) {
  
  accepted.formats.emat <- c("matrix", "Matrix", "data.frame")  

  # ensure datamatrix is a matrix
  if (!is.matrix(datamatrix)) {
    if (class(datamatrix) %in% accepted.formats.emat) {
      datamatrix <- as.matrix(datamatrix)
    } else {
      stop(paste("datamatrix needs to be in one of the following formats:", paste(accepted.formats.emat, collapse = "; ")))
    }    
  }
  if (is.null(colnames(datamatrix))) { colnames(datamatrix) <- as.character(1:ncol(datamatrix)) }
  if (is.null(rownames(datamatrix))) { rownames(datamatrix) <- as.character(1:nrow(datamatrix)) }  

  datamatrix
}



#' filter.network
#' 
#' 
#' Include to the network only the edges with the highest mutual information,
#' calculated based on the first principal components.
#' 
#' @param network network
#' @param delta associated cost function value changes for each node merge
#' @param datamatrix datamatrix
#' @param params parameters
#' @return Filtered network
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @examples
#' 
#' # TBA
#' 
filter.network <- function (network, delta, datamatrix, params) {

  # Include at maximum speedup.max.edges for each node in the network
  # based on mutual information
  # return list of edges to keep

  remove.edges <- c()
  uniq.edges <- unique(network[1,])  

  for (idx in 1:length(uniq.edges)) {
    if (params$verbose) {message(paste(idx, "/",  length(uniq.edges)))}
    a <- uniq.edges[[idx]]

    # Pick edges for this node
    eds <- which(network[1, ] == a)

    # Calculate MI scores
    #require(minet)
    mis <- c()
    mi.cnt <- 0  
    # FIXME: could be parallelized
    for (edge in eds){ # which(is.na(delta))
      mi.cnt <- mi.cnt + 1

      # Pick node indices
      i <- network[2, edge]

      # For singletons
      dat <- cbind(datamatrix[, a], datamatrix[, i])
                   
      mis[[mi.cnt]] <- build.mim(dat, estimator="mi.empirical", disc = "equalwidth", nbins = params$nbins)[1, 2]
    }
    
    # Edges to remove
    remove.edges <- eds[na.omit(order(mis, decreasing = TRUE)[-seq(params$speedup.max.edges)])]
    keep.edges <- setdiff(1:ncol(network), remove.edges)

    # Filter out remove.edges
    if (length(remove.edges) > 0) {
      network <- network[, keep.edges]
      delta <- delta[keep.edges] 
    }

  }

  list(network = network, delta = delta)

}

edge.delta <- function (edge, network, network.nodes, datamatrix, params, model.nodes, model) {

    if ( params$verbose ) { message(paste('Computing delta values for edge ', edge, '/', ncol(network), '\n')) }
    a <- network[1, edge]
    b <- network[2, edge]
    vars <- network.nodes[c(a, b)]

    tmp <- mixture.model(matrix(datamatrix[, vars], nrow( datamatrix )), vars, params) 
    model <- tmp$model # FIXME: perhaps the 'model' is not needed when model.params is given. Check and remove.
    model.params <- tmp$params

    # Compute COST-value for two independent subnets vs. joint model
    # Negative free energy (-cost) is (variational) lower bound for P(D|H)
    # Use it as an approximation for P(D|H)
    # Cost for the indpendent and joint models
    # -cost is sum of two independent models (cost: appr. log-likelihoods)
    costind.ab     <-  info.criterion(model.nodes[[a]]$Nparams + model.nodes[[b]]$Nparams, params$Nlog, -(model.nodes[[a]]$free.energy + model.nodes[[b]]$free.energy), criterion = params$information.criterion)
    costjoint.ab   <-  info.criterion(model.params$Nparams, params$Nlog, -model.params$free.energy, criterion = params$information.criterion)

    # NOTE: COST is additive so summing is ok
    # change (increase) of the total COST / cost
    delt <- as.numeric(costjoint.ab - costind.ab)

    # Store these only if it would improve the cost; otherwise never needed
    if (-delt > params$merging.threshold) {
      mod.pair <- model.params
    } else {
      mod.pair <- 0
    }
			
    return(list(c(mod.pair, delt = delt)))
}



build.mim <- function (dataset, estimator = "spearman", disc = "none", nbins = sqrt(NROW(dataset))) 
{
    # This function is licensed under cc-by-sa 3.0
    # Modified from minet 3.6.0 to use discretize correctly (not exported in build.mim which may cause occasional function name conflicts)

    if (disc == "equalfreq" || disc == "equalwidth" || disc == 
        "globalequalwidth") 
        dataset <- discretize(dataset, disc, nbins)
    if (estimator == "pearson" || estimator == "spearman" || 
        estimator == "kendall") {
        mim <- cor(dataset, method = estimator, use = "complete.obs")^2
        diag(mim) <- 0
        maxi <- 0.999999
        mim[which(mim > maxi)] <- maxi
        mim <- -0.5 * log(1 - mim)
    }
    else if (estimator == "mi.mm") 
        estimator = "mm"
    else if (estimator == "mi.empirical") 
        estimator = "emp"
    else if (estimator == "mi.sg") 
        estimator = "sg"
    else if (estimator == "mi.shrink") 
        estimator = "shrink"
    else stop("unknown estimator")
    if (estimator == "mm" || estimator == "emp" || estimator == 
        "sg" || estimator == "shrink") {
        mim <- mutinformation(dataset, method = estimator)
        diag(mim) <- 0
    }
    mim[mim < 0] <- 0
    mim
}

#build.mim.c <- cmpfun(build.mim)

#######################################

discretize  <- function (X, disc = "equalfreq", nbins = sqrt(NROW(X))) 
{
    # This function is licensed under cc-by-sa 3.0
    # Modified from minet 3.6.0 internal function

    X <- as.data.frame(X)
    varnames <- names(X)
    dimensions <- dim(X)
    X <- data.matrix(X)
    dim(X) <- dimensions
    res <- NULL
    if (disc == "equalfreq") 
        res <- .Call("discEF", X, NROW(X), NCOL(X), as.integer(nbins), 
            DUP = FALSE, PACKAGE = "infotheo")
    else if (disc == "equalwidth") 
        res <- .Call("discEW", X, NROW(X), NCOL(X), as.integer(nbins), 
            DUP = FALSE, PACKAGE = "infotheo")
    else if (disc == "globalequalwidth") 
        res <- as.vector(cut(X, nbins, labels = FALSE))
    else stop("unknown discretization method")
    dim(res) <- dimensions
    res <- as.data.frame(res)
    names(res) <- varnames
    res
}





#' pick.model.parameters
#' 
#' @param m vdp.mixt output
#' @param nodes node names for naming purposes 
#' @return Model parameters
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @export
#' @examples #
pick.model.parameters <- function (m, nodes) {
 
  # Pick parameters
  w    <- m$posterior$weights    # component weights
  mu   <- m$posterior$centroids  # component centroids
  sds  <- m$posterior$sds        # component standard devs

  rownames(mu) <- rownames(sds) <- names(w) <- paste("Mode", 1:length(w), sep = "-")
  colnames(mu) <- colnames(sds) <- nodes		       

  # For mu and std, rows correspond to the mixture components, in w the elements
  list(mu = mu, sd = sds, w = w, free.energy = m$free.energy, Nparams = m$posterior$Nparams)

}


get.subnet <- function (res, subnet.id) {

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", subnet.id, sep="")
  }
    
  # Nodes for a given subnet
  get.subnets(res)[[subnet.id]]

}


compute.weight <- function (pt, mu, vars, xt) {

  # pt <- qOFz[t,] # P(c|t) for all c
  # xt  <- dat[t, ] # data point
  
  # Initial weights are zero
  w <- rep.int(0, nrow(mu))  
  sds <- sqrt(vars)

  # Avoid overflows: no probability can be 0 exactly. Add negligible constant in these cases
  pt[pt < 1e-320] <- 1e-320
  pt <- pt/sum(pt) # renormalize

  # set arbitrary weigh for the first cluster
  # (only relations between weights matter)
  # normalize later
  #w[[1]] <- 1 # added below
  # operate on log domain to avoid floating errors

  logdens <- sum(dnorm(xt, mu[1, ], sds[1, ], log = TRUE))
 
  if (nrow(mu) > 1) {
    logw <- sapply(2:nrow(mu), function (i) { 
      #print(prod(dnorm(xt, mu[i, ], sds[i, ])))
  
      #dens * pt[[i]]/(pt[[1]] * prod(dnorm(xt, mu[i, ], sds[i, ]))) 

      logdens + log(pt[[i]]) - log(pt[[1]]) - sum(dnorm(xt, mu[i, ], sds[i, ], log = TRUE))

    })
    w <- exp(logw)

    # Densities in original domain, standardized s.t. first weight is 1
    w <- c(1, w)

    # normalize to unity and return
    return(w/sum(w))
  } else {
    return(1) #w <- 1
  }
    
}



############################################


retrieve.model <- function (model, subnet.id) {
  # Note: the algorithm has some stochasticity in initialization etc.
  # so the results may not be exactly same each time; but they should
  # be sufficiently similar in any case
  
  # Former: get.model

  # model: output from run.netresponse function
  # subnet.id: id/index of the subnet to check
  # level: which agglomeration step
  # datamatrix for which the model was calculated
  
  #  Copyright (C) 2008-2011 Leo Lahti
  #  Licence: GPL >=2

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", subnet.id, sep="")
  }
  
  # Get subnet nodes
  nodes <- model@subnets[[subnet.id]]

  # Compute the model
  vdp.mixt(matrix(model@datamatrix[, nodes], nrow(model@datamatrix)))

}



# INPUT:   data, hp_posterior, hp_prior, opts
# OUTPUT:  list(free_energy,hp_posterior,data,c)
#        * free_energy: free energy of the best split found
#        * hp_posterior: posterior info of the best split found
#        * c: index of the cluster that resulted in the best split found
#
# DESCRIPTION: Implements the VDP algorithm steps 2 to 4.

find.best.splitting <- function(data, hp.posterior, hp.prior, opts, min.size = 5){

  # min.size: minimum size of a component required for splitting
  
  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner  

  dat <- data$given.data$X1
  epsilon <- 1e-10    # FIXME: should this be a tunable function parameter?
  c.max <- opts$c.max 

  # ALGORITHM STEP 2
 
  # Sort clusters by size, use at most c.max candidates and ensure cluster size is > 2
  candidates <- order(hp.posterior$Nc, decreasing = TRUE)
  candidates <- candidates[hp.posterior$Nc[candidates] > 2]
  # candidates <- which(hp.posterior$Nc > 2)
  if ( length(candidates) == 0 ) { c <- 1 }

  qOFz <- mk.qOFz(data, hp.posterior, hp.prior, opts)
  fc   <- mk.E.log.q.p.eta(data, hp.posterior, hp.prior, opts)
  log.lambda <- mk.log.lambda(data, hp.posterior, hp.prior, opts) 
  sub.data  <- data

  # ALGORITHM STEP 3 (3a,3b,3c) check which split gives best improvements in cost

  new.qOFz.list <- list()   #Initialize
  new.free.energy <- rep(Inf, min(c.max, length(candidates)))

  for (c in candidates[1:min(c.max, length(candidates))]) {

    # relating.n has the indexes of data points that belonged to the candidate cluster prior
    # to splitting (that's why it is the sum over the now 2 clusters (after splitting).
    # REMARK: Is 0.5 ok? - when there are lots of clusters it is natural to 
    # assume some points will have less than 0.5 for any cluster.

    relating.n <- which(qOFz[, c] > 0.5)
    if (length(relating.n) == 0) { next } else {}

    # ALGORITHM STEP 3a. split the candidate cluster
    # Split cluster c in qOFz into two smaller ones

    new.c     <- ncol(qOFz)
    new.qOFz  <- split.qofz(qOFz, c, new.c, dat, opts$speedup, min.size)

    new.K     <- ncol( new.qOFz )
    sub.qOFz  <- new.qOFz[relating.n, unique(c(c, new.c, new.K))]
    
    # Ensure it remains a matrix
    sub.data$given.data$data <- sub.data$given.data$X1 <- array(dat[relating.n, ],
                                    dim = c(length(relating.n), ncol(dat)))    

    # ALGORITHM STEP 3b and 3c
    # update the posterior of the split clusters for a small number of iter.
    # update_posterior sorts clusters by size

    sub.hp.posterior <- mk.hp.posterior(sub.data, sub.qOFz, hp.prior, opts)
    dummylist        <- updatePosterior(sub.data, sub.hp.posterior, hp.prior, opts, 10, 0)
    sub.hp.posterior <- dummylist$hp.posterior
    sub.qOFz         <- dummylist$qOFz

    # FIXME: check this already previously for c == new.c? 
    if(ncol( sub.qOFz ) < 3) { next } else { } 

    # If there are more than 1 empty components then go to next step
    if(sum(colSums(sub.qOFz) < epsilon) > 1) { next } else { }

    sub.log.lambda <- mk.log.lambda(data, sub.hp.posterior, hp.prior, opts)
    insert.indices <- c(c, new.c, new.K:(new.K + ncol(sub.qOFz) - 3))
    
    if(max(insert.indices) > ncol(log.lambda)){
      new.log.lambda <- cbind(log.lambda, array(0, dim = c(nrow(log.lambda), max(insert.indices) - ncol(log.lambda))))
    } else{ new.log.lambda <- log.lambda }

    
    new.log.lambda[, insert.indices] <- sub.log.lambda
    new.fc <- fc
    new.fc[insert.indices] <- mk.E.log.q.p.eta(sub.data, sub.hp.posterior, hp.prior, opts)
    new.free.energy[[c]] <- mk.free.energy(data, sub.hp.posterior, hp.prior, opts, new.fc, new.log.lambda)$free.energy

    # if new.qOFz is not large enough to accommodate update from sub.qOFz then add columns
    if (ncol(new.qOFz) < max(insert.indices)) {
      new.qOFz <- cbind(new.qOFz, array(0, dim = c(nrow(new.qOFz), max(insert.indices) - ncol(new.qOFz))))
    }
    
    new.qOFz[relating.n, ] <- 0
    new.qOFz[relating.n, insert.indices] <- sub.qOFz
    new.qOFz.list[[c]]     <- new.qOFz
  }

  # Select cluster split that minimizes free energy
  c <- which.min(new.free.energy)
  free.energy <- new.free.energy[[c]]  

  
  if(is.infinite(free.energy)){
    c <- -1
  } else { hp.posterior <- mk.hp.posterior(data, new.qOFz.list[[c]], hp.prior, opts) }

  list( free.energy = free.energy,
       hp.posterior = hp.posterior,
               data = data,
                  c = c)
}



#find.best.splitting.c <- cmpfun(find.best.splitting)

########################################



find.best.neighbor <- function (G, max.subnet.size, network, delta) {

  # Order edges by delta values. The two subnets with the smallest delta are
  # joined unless the merged subnet exceeds max size.
  o <- order(delta)

  # Check size of the resulting merged subnet one-by-one, starting from the smallest
  best.found <- FALSE
  cnt <- 0
  best.edge <- a <- b <- NULL
  mindelta <- Inf
  while (!best.found) {
    cnt <- cnt + 1
    ind <- o[[cnt]]
    z <- network[1, ind]
    i <- network[2, ind]    

    # Finish when new merged subnetwork that does not exceed max size is found
    if (length(c(G[[z]], G[[i]])) <= max.subnet.size){
      # Sort a and b
      a <- min(c(z, i))
      b <- max(c(z, i))
      best.edge <- ind
      mindelta <- delta[[best.edge]]
      best.found <- TRUE
    }  else {
      # Put infinite cost for merges that would exceed maximum size
      delta[[ind]] <- Inf 
    }
  }

  list(a = a, b = b, mindelta = mindelta, best.edge = best.edge, delta = delta)

}


#find.best.neighbor.c <- cmpfun(find.best.neighbor)



join.subnets <- function (network, delta, best.edge) {
  # for edge matrices
  a <- network[1, best.edge]
  b <- network[2, best.edge]  
  
  # replace b nodes by a: this in effect transfers b edges to a edges
  # and removes b completely
  network[network == b] <- a
  
  # remove delta values for nodes associated with a node
  # since this node now has different (merged) set of features
  # and model needs to be recalculated
  inds <- apply(network == a, 2, any)
  delta[inds] <- NA

  # sort entries (row1 < row2)
  network <- apply(network, 2, sort)
      
  list(network = matrix(network, 2), delta = delta)
}

#join.subnets.c <- cmpfun(join.subnets)

########################################################

check.bins <- function (difexp, mybreaks) {

  # check color scale bin for each expression value
  bins <- c()
  for (i in 1:length(difexp)) {
      # which color bins are smaller than our difexp value
      # (for probet: i, mode:mode)
      inds <- which(difexp[[i]] > mybreaks)
      if (length(inds) == 0) {
        bins[[i]] <- 1
      } else if (length(inds) > 0)  {
        bins[[i]] <- max(inds) + 1
      }
  }

  bins
}



###############################################################################


# INPUT: matrix (q_of_z)
# OUTPUT: matrix

# DESCRIPTION: Sorted matrix in decreasing fashion based on the value
#              of colSums.  Remark: The last column of the matrix is
#              kept in place (it is not sorted).

sortqofz <- function(qOFz){

  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner
  
  Nc <- colSums(qOFz)
  I  <- order(Nc[-length(Nc)], decreasing = TRUE) 
  I  <- c(I, ncol(qOFz)) # add last column element separately (outside of ordering)
                                    
  # Order the cols and ensure qOFz remains a matrix
  array(qOFz[, I], dim = dim(qOFz))
  
}

#sortqofz.c <- cmpfun(sortqofz)

############################################################


# INPUT:   data   - matrix with data vectors
#          K      - number of clusters
# OUTPUT:  q_of_z - matrix of size N*(K+1).
# DESCRIPTION: This function assigns data randomly to K clusters by drawing cluster
#              membership values from a uniform distribution.

rand.qOFz <- function(N, K){

  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  qOFz <- matrix( runif(N*(K+1)), N, K + 1)
  qOFz[, K + 1] <- 0

  # normalize and return (each row should sum up to 1)
  qOFz/rowSums(qOFz) 

}

#rand.qOFz.c <- cmpfun(rand.qOFz)

############################################################


# INPUT: "old" free_energy value, "new" free_energy value, options
# OUTPUT: bool: 0 if the there was no significant improvement.
#               1 if new_free_energy is smaller than free_energy (more than opts$threshold).



free.energy.improved <- function(free.energy, new.free.energy,
                                 warn.when.increasing, threshold)
{

  #INPUT: "old" free.energy value, "new" free.energy value, options
  #OUTPUT: bool: 0 if the there was no significant improvement.
  #              1 if new.free.energy is smaller than free.energy
  #                (more than opts$threshold).

  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  diff <- new.free.energy - free.energy

  v <- abs(diff/free.energy)
  
  if(is.nan(v) || v < threshold){
    bool <- 0
  } else { 
    if(diff > 0){
      if( warn.when.increasing ){
        if( v > 1e-3 ){
          stop(c("the free energy increased. The diff is ", toString(diff)))
        } else {
          warning(c("the free energy increased. The diff is ", toString(diff)))
        }
      }
      bool <- 0
    } else {
      bool <- ifelse(diff == 0, 0, 1)
    }
  }
  
  bool
}

#################################################################################


# INPUT:   data: structure with data matrix
#          hp_posterior: posterior information for the current mixture model
#          hp_prior: prior information
#          opts: options list.
#
# OUTPUT:  list(free_energy, hp_posterior, hp_prior, data);
# 
# DESCRIPTION: Read the main description on the beginning of the file.


greedy <- function(data, hp.posterior, hp.prior, opts, min.size){

  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner
  
  free.energy <- mk.free.energy(data, hp.posterior, hp.prior, opts)$free.energy

  while(1){

    # ALGORITHM STEP 2-4

    templist <- find.best.splitting(data, hp.posterior, hp.prior, opts, min.size)
    new.free.energy  <- templist$free.energy
    new.hp.posterior <- templist$hp.posterior
    c                <- templist$c
    if ( c == (-1) ) { break } # infinite free energy -> break splitting

    # ALGORITHM STEP 5
    dummylist <- updatePosterior(data, new.hp.posterior, hp.prior,
                                  opts, ite = opts$ite, do.sort = 1)

    new.free.energy  <- dummylist$free.energy
    new.hp.posterior <- dummylist$hp.posterior
    
    # ALGORITHM STEP 6
    if( free.energy.improved(free.energy, new.free.energy, 0, opts$threshold) == 0 ) {
      break #free.energy didn't improve, greedy search is over
    } 

    free.energy  <- new.free.energy
    hp.posterior <- new.hp.posterior

  }

  list(free.energy = free.energy,
       hp.posterior = hp.posterior,
       hp.prior = hp.prior,
       data = data)
}

###############################################################################


#INPUT:   data, hp_posterior, hp_prior, opts, ite, do_sort
#             * do_sort: TRUE/FALSE: indicates whether the clusters should be 
#                        sorted by size or not.
#             * ite: number of update iterations.
#OUTPUT:  Updated parameters: list(free_energy, hp_posterior, q_of_z)
#DESCRIPTION: Updates the posterior of the mixture model. if do_sort=true it also
#             sorts the cluster labels by size. (i.e. cluster 1 = largest cluster)


updatePosterior <- function(data, hp.posterior, hp.prior, opts, ite = Inf, do.sort = 1) {
                            
  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  epsilon <- 1e-10
  free.energy <- Inf
  i <- last.Nc <- start.sort <- 0
  opts.internal <- opts
  break.loop <- FALSE
  
  while( 1 ){ 
    i <- i + 1

    new.free.energy <- Inf
    cnt <- 0
    while (is.infinite(new.free.energy) && cnt <= 10) {
      
      templist <- mk.free.energy(data, hp.posterior, hp.prior, opts.internal)
      new.free.energy <- templist$free.energy
      log.lambda <- templist$log.lambda
      if( is.infinite(new.free.energy) ) {
        warning("Free energy not finite: adding implicit noise.")
        opts.internal$implicitnoisevar <- opts.internal$implicitnoisevar + 0.1
      }
      cnt <- cnt + 1
    }

    if ( (is.finite(ite) && i>= ite) ||
         (is.infinite(ite) && 
          free.energy.improved(free.energy, new.free.energy, 0, opts.internal$threshold) == 0)){
        free.energy <- new.free.energy
        if( do.sort && opts.internal$do.sort && (!start.sort) && is.finite(free.energy)){
          start.sort <- 1
        } else break # this will break the while loop
      }

      last.Nc <- hp.posterior$Nc
      free.energy <- new.free.energy
      qOFz <- mk.qOFz(data, hp.posterior, hp.prior, opts.internal, log.lambda)

      # if the last component is not 'empty' and max number components not reached
      # add a new empty component
      if(sum(qOFz[, ncol(qOFz)]) >= epsilon && ncol(qOFz) < opts$c.max){ 
        qOFz <- cbind(qOFz, 0) 
      }

      # Sort components by size (note: last component kept in its place)
      if( start.sort ){ qOFz <- sortqofz(qOFz) }

      # Pick at most c.max+1 components, remove the smallest one, except the one added in this iteration (ie. c.max + 1)
      # FIXME: consider how to improve the implementation regarding to this!
      # (3/2012)
      #if (ncol(qOFz) == opts$c.max + 1) {
      # 	 qOFz <- qOFz[, setdiff(1:ncol(qOFz), opts$c.max), drop = FALSE]        
      #} 

      # If the smallest of the previous components
      # (excluding the one added in this interation)
      # is empty then remove it
      # (if new component was not added this is ok to have as well)
      if(sum(qOFz[, ncol(qOFz) - 1 ]) < epsilon){
        qOFz <- matrix(qOFz[, -(ncol(qOFz) - 1)], nrow(qOFz))
      }

      qOFz <- matrix(qOFz/rowSums(qOFz), nrow(qOFz)) # probabilities sum to one      
      hp.posterior <- mk.hp.posterior(data, qOFz, hp.prior, opts.internal)    

  }

  list( free.energy = free.energy,
       hp.posterior = hp.posterior,
               qOFz = qOFz)
}


sumlogsumexp <- function(log.lambda){.Call("vdpSumlogsumexp", log.lambda, PACKAGE = "netresponse")}


softmax <- function( A ){
  qOFz <- .Call("vdpSoftmax", A, PACKAGE = "netresponse")

  # Ensure that this remains a matrix even if it has 1-dimension on rows or cols
  qOFz <- array(qOFz, dim = dim(A))
  
  # In rare (< 1 / 1e6?) cases, and particularly when sample size is
  # large (>1500) the vdpSoftmax function seems to produce NaNs. This
  # occurred for particularly low values of A: A[,1] in the range <
  # -800 and below while other values in A[,1] were at least -400 and
  # mostly in range -100 ... 0; for A[,2] typical range around -800,
  # and for the NaN outlier: -1000. In both cases it was A[209,1] and
  # A[209,2] i.e. the same sample.  This is so rare that the effect on
  # the results is negligible, but this should be diagnosed and fixed
  # asap. As this is part of iteration, simply replace the defected sample with
  # equal probability in all groups.
  if (sum(!is.na(qOFz)) > 0) {
    # FIXME: optimize with apply!
    # Detect empty components and ignore
    inds <- c()
    for (i in 1:ncol(qOFz)) {
      if (sum(na.omit(qOFz[, i])) == 0) {
        inds <- c(inds, i)
        qOFz[, i] <- 0
      }
    }

    for (i in 1:nrow(qOFz)) {
      if (sum(is.na(qOFz[i, ])) > 0) {
        inds2 <- setdiff(seq(ncol(qOFz)), inds)
        qOFz[i, inds2] <- rep.int(1/length(inds2), length(inds2))
      }
    }
  }
  
  qOFz
  
}

mk.hp.posterior <- function(data, qOFz, hp.prior, opts){

  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  dat <- data$given.data$X1
  
  # Ensure that qOFz is a matrix
  qOFz <- matrix(qOFz, nrow(dat))

  # If qOFz exceeds max cluster size then remove one cluster 
  # (the second last one, which assumes clusters are sorted and the last is a new one)
  # FIXME 3/2012 quick hack - consider better implementations regarding this
  if (ncol(qOFz) > opts$c.max + 1) {
    inds <- setdiff(1:ncol(qOFz), ncol(qOFz) - 1)
    qOFz <- matrix(qOFz[, inds], nrow(dat))
    qOFz <- matrix(qOFz/rowSums(qOFz), nrow(dat))
  }

  # Compatibility variables not needed for the current functionality
  tmp.realS <- X2 <- dimX2 <- 0

  out <- .Call("mHPpost",
               dat,
               ncol(dat),
               nrow(dat),
               X2, dimX2,
               tmp.realS, opts$implicitnoisevar,
               hp.prior$Mumu, hp.prior$S2mu,
               hp.prior$AlphaKsi, hp.prior$BetaKsi,
               hp.prior$U.p, hp.prior$alpha,
               qOFz, ncol(qOFz), PACKAGE = "netresponse")

  qOFz <- matrix(out$qOFz, nrow(qOFz))

  #if (ncol(qOFz) > opts$c.max) {
  #}

  hp.posterior <- list(
    Mubar     = matrix(out$Mubar,    ncol(qOFz)),
    Mutilde   = matrix(out$Mutilde,  ncol(qOFz)),
    KsiAlpha  = matrix(out$KsiAlpha, ncol(qOFz)),
    KsiBeta   = matrix(out$KsiBeta,  ncol(qOFz)),
    gamma     = matrix(out$gamma, 2),
    Nc        = out$Nc,
    qOFz      = qOFz,
    Uhat      = out$Uhat)

  hp.posterior
}


mk.hp.prior <- function(data, opts){

  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

 # INPUT:   data     - matrix with data vectors
 #          opts     - list with algorithm options
 # OUTPUT:  hp_prior - list with prior information
 #
 # DESCRIPTION: Prior for the mean     = mean of data
 #              Prior for the variance = variance of data


  dat <- data$given.data$X1 # real-valued. Data to be clustered.

  Mean  <- colMeans(dat)     # mean of each dimension
  Var   <- colVariances(dat, Mean) # Variance of each dimension
                  #colSums((dat - rep(Mean, each = nrow(dat)))^2)/nrow(dat) 

  # priors for distribution of codebook vectors Mu ~ N(MuMu, S2.Mu)..
  #list(Mumu = Mean, S2mu = Var, U.p = Inf)
  # priors for data variance Ksi ~ invgam(AlphaKsi, BetaKsi)
  # variance is modeled with inverse Gamma distribution
  # FIXME: some of these are redundant, remove to save memory
  list(Mumu = Mean, S2mu = Var, U.p = Inf, AlphaKsi = rep(opts$prior.alphaKsi, ncol(dat)), BetaKsi = rep(opts$prior.betaKsi, ncol(dat)), alpha = opts$prior.alpha)
  
}


########################################################################################


# INPUT:   data, hp_posterior, hp_prior, opts
# OUTPUT:  free_energy: value of mixture model's free energy
#          log_lambda: Used for posterior of labels q_of_z <- softmax(log_lambda);
# DESCRIPTION: ...


mk.free.energy <- function(data, hp.posterior, hp.prior, opts,
                           fc = NULL, log.lambda = NULL)
{

  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner
  
  if( is.null(fc) || is.null(log.lambda) ){
    fc <- mk.E.log.q.p.eta(data, hp.posterior, hp.prior, opts)      # 1*K
    log.lambda <- mk.log.lambda(data, hp.posterior, hp.prior, opts) # N*K
  } 

  hpgsum <- colSums(hp.posterior$gamma)
  dig <- digamma(hpgsum)
  
  E.log.p.of.V <- lgamma(hpgsum) -
      lgamma(1 + hp.prior$alpha) -
      colSums(lgamma(hp.posterior$gamma)) +
      lgamma(hp.prior$alpha) +
      ( (hp.posterior$gamma[1, ] - 1) * 
        (digamma(hp.posterior$gamma[1, ]) - dig)) +
      ( (hp.posterior$gamma[2, ] - hp.prior$alpha) *
        (digamma(hp.posterior$gamma[2, ]) - dig))

  extra.term <- sum(E.log.p.of.V)
  free.energy <- extra.term + sum(fc) - sumlogsumexp(log.lambda)

  # Return
  list(free.energy = free.energy, log.lambda = log.lambda)

}
  
#mk.free.energy.c <- cmpfun(mk.free.energy)

mk.qOFz <- function(data, hp.posterior, hp.prior, opts, log.lambda = NULL){

  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  if( is.null(log.lambda) ){
    log.lambda <- mk.log.lambda(data, hp.posterior, hp.prior, opts)
  }

  qOFz <- softmax( log.lambda )

  # Do not allow empty clusters
  as.matrix(qOFz[, !colSums(qOFz) == 0], nrow(log.lambda))
  
}



mk.log.lambda <- function(data, hp.posterior, hp.prior, opts){


  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  dat <- data$given.data$X1
  
  # Ensure matrix for data
  if (!is.matrix(dat)) {stop("Error in mk.log.lambda: dat should be a matrix!")}
  
  # Compatibility variables, not needed for current functionality
  tmp.realS <- X2 <- dimX2 <- 0

  out <- .Call("mLogLambda",
               dat,
               ncol(dat),
               nrow(dat),
               X2, dimX2,
               tmp.realS, opts$implicitnoisevar,
               hp.prior, hp.posterior,
               PACKAGE = "netresponse")
  
  matrix(out, nrow(dat))
  
}


###########################################################################


# INPUT:   data, hp_posterior, hp_prior, opts
# OUTPUT:  matrix [1xk]: used to compute the free_energy formula. 
# DESCRIPTION: Regards the gaussian model's parameters.



mk.E.log.q.p.eta <- function(data, hp.posterior, hp.prior, opts){

  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner
  
  # returns E [ log q(eta)/p(eta) ].q
  # fc: 1 by k

  dat <- data$given.data$X1
  
  # Ensure matrix for data
  if (!is.matrix(dat)) {stop("Error in mk.E.log.q.p.eta: dat is not a matrix!")}
  
  N  <- nrow(dat)
  M1 <- ncol(dat)
  K  <- nrow(hp.posterior$Mubar)
  l.codebook <- (- M1/2) * matrix(1, K)
  Ksi.log <- (digamma(hp.posterior$KsiAlpha) - log(hp.posterior$KsiBeta))
  
  for (j in 1:M1) {
    l.codebook <- l.codebook +
      .5 * ( log(hp.prior$S2mu[[j]] / hp.posterior$Mutilde[, j]) +
             ( (hp.posterior$Mubar[, j] - hp.prior$Mumu[[j]])^2 +
                hp.posterior$Mutilde[, j] 
             ) / hp.prior$S2mu[[j]]
           ) +
      lgamma(hp.prior$AlphaKsi[[j]]) - 
      lgamma(hp.posterior$KsiAlpha[, j]) +
      hp.posterior$KsiAlpha[, j] * log(hp.posterior$KsiBeta[, j]) -
      hp.prior$AlphaKsi[[j]] * log(hp.prior$BetaKsi[[j]]) + 
      (hp.posterior$KsiAlpha[, j] - hp.prior$AlphaKsi[[j]]) * Ksi.log[, j] +
      (hp.prior$BetaKsi[[j]] - hp.posterior$KsiBeta[, j]) *
        (hp.posterior$KsiAlpha[, j] / hp.posterior$KsiBeta[, j])
  }

  t(l.codebook)
}

#mk.E.log.q.p.eta.c <- cmpfun(mk.E.log.q.p.eta)

#################################################


colVariances <- function (dat, Mean) {
  # This is about 5x faster than apply(dat, 2, var)
  #max(abs(colVariances(dat) - apply(dat,2,var)))
  colSums((dat - rep(Mean, each = nrow(dat)))^2)/(nrow(dat) - 1)
}

#colVariances.c <- cmpfun(colVariances)

############################################################################

# INPUT:   data, qOFz, hp_posterior, hp_prior, opts
# OUTPUT:  list(new.qOFz, new.c);
#             * new.qOFz: posterior over labels including the split clusters.
#             * new.c: index of the newly created cluster.
# DESCRIPTION: Implements the VDP algorithm step 3a.

split.qofz <- function(qOFz, c, new.c, dat, speedup = TRUE, min.size = 4){

  # compute the first principal component of the candidate cluster,
  # not the whole data set.

  # min.size option sets the required minimum size of a cluster for
  # splitting; smaller clusters are not splitted

  # Pick sample indices and samples corresponding to cluster c
  cluster_assignments <- apply(qOFz, 1, which.max);
  indices <- which(cluster_assignments == c);
  if (length(indices) < min.size) {
    #"Component must have at least min.size samples to be splitted."
    # -> no splitting
    new.qOFz <- qOFz
  } else {

    component.data <- matrix(dat[indices,], length(indices))
  
    # If the number of samples is high calculating PCA might take long
    # but can be approximated by using less samples:

    pcadata <- component.data
    
    if ( speedup ) {

      # when a candidate cluster, C, is split to generate two new
      # clusters, it is split by mapping the data onto the first
      # principal component of the data in C and then splitting that in
      # half. To speed up, one can compute an approximate first
      # principal component by considering a randomly selected subset of
      # the data belonging to C, and computing its first principal
      # component.

      # number of samples in this component
      ns <- nrow(component.data) 
      nd <- ncol(component.data) 

      # If component size exceeds cmax, 
      # use only a random subset of data to calculate PCA
      # size of the random subset increases slowly (linearly) 
      # with component size. 
      cmax <- 20 #take at least this many random samples    
      nr <- min(ns, cmax + floor(sqrt(ns))) # do not take more samples than are available
      rinds <- sample(ns, nr)

      # Pick random subset of the component data and accompanying indices
      # to speed up PCA calculations
      pcadata <- matrix(component.data[rinds,], nrow = nr)
      indices <- indices[rinds]
    }

    # Split the cluster based on the first PCA component
    # FIXME: compare speed with other PCA implementations and select fastest
    dir <- prcomp(pcadata)$x[,1]
    I1 <- indices[dir >= 0];
    I2 <- indices[dir < 0];

    # Initialize split by adding a zero column on qOFz

    # If one of qOFz clusters is empty, then do not create new clusters but instead fill in the empty cluster
    # during cluster split.
    # FIXME: ensure already in creating qOFz-matrices that no zero columns are allowed. This will
    # avoid the need to address the issue here.
    # -> OK, done this. w remove this unnecessary check here and test if
    # the code works ok
    empty.cols <- (colSums(qOFz) == 0)
    if ( !any(empty.cols) ) { # no empty columns -> add an empty cluster
      new.qOFz <- array(0, dim = c(nrow(qOFz), ncol(qOFz) + 1))
      new.qOFz[,  -new.c] <- qOFz
    } else { # an empty column -> no need to add new clusters
      new.qOFz <- qOFz
      new.c <- which(empty.cols)[[1]]
    }

    # Split this component (samples given in I1, I2) into two smaller components
    new.qOFz[ I1, c]     <- qOFz[ I1, c]
    new.qOFz[ I2, c]     <- 0 # Remove entries from cluster c
    new.qOFz[ I2, new.c] <- qOFz[ I2, c] # Add same entries to cluster new.c
  }

  new.qOFz

}

#split.qofz.c <- cmpfun(split.qofz)

###############################################################################


#fsort <- function (df, sortvar, decreasing = FALSE) {
fsort <- function (df, sortvar) {

  o <- order(df[[sortvar]])
  #if (decreasing) {o <- rev(o)}
  df[o,]

}

#fsort.c <- cmpfun(fsort)

##################################################################
	
plotMatrix.2way <- function (mat, mybreaks = NULL, maintext = "", xlab = "", ylab = "", mypalette = NULL, interval = .1, cex.main = 1, xaxis = FALSE, yaxis = TRUE, row.tick = 1, col.tick = 1, cex.xlab = .9, cex.ylab = .9, cex.lab = .9, limit.trunc = 0, mar = c(5, 4, 4, 2), ...) {

  # mat: differential expression matrix to plot in two-color palette
  # interval: interval for palette color switches
  # FIXME: synchronize with PlotMatrix in sorvi package  
    
  require(graph)
  #require(RBGL)
  require(Rgraphviz)
  require(graphics)
	   
  if (length(mybreaks) == 0)  {
    m <- max(round(max(abs(mat)), limit.trunc) - interval, 0)
    mm <- m + interval/2
    vals <- seq(interval/2,mm,interval)
    # Set breaks evenly around zero
    mybreaks  <- c(-(m+1e6),c(-rev(vals),vals),m+1e6)
  }
		  
  if (length(mypalette)==0) {
    mypalette <- colorRampPalette(c("blue", "black", "red"),space = "rgb")
    my.colors <- mypalette(length(mybreaks)-1)
  } else {
    my.colors <- mypalette(length(mybreaks)-1)
  }
		      
  # transpose and revert row order to plot matrix in the same way it
  # appears in its numeric form
  par(mar = mar)
  image(t(mat[rev(seq(nrow(mat))),]), col = my.colors, xaxt='n', yaxt='n', zlim=range(mybreaks), breaks=mybreaks, main=maintext, xlab=xlab, ylab=ylab, cex.lab = cex.lab, cex.main = cex.main)

  if (yaxis) {
      
    v <- seq(1, nrow(mat), row.tick) # take every nth index
    axis(2, at = seq(0,1,length = nrow(mat))[v], labels = rev(rownames(mat))[v], 
    	    cex.axis=cex.ylab, las = 2)
    
  }
  
  if (xaxis) {    

    v <- seq(1, ncol(mat), col.tick) # take every nth index
    axis(1, at = seq(0,1,length = ncol(mat))[v], labels = colnames(mat)[v], 
    	    cex.axis = cex.xlab, las=2)

  }
    
  return(list(palette = my.colors, breaks = mybreaks))
      	  
}


check.bins <- function (difexp, mybreaks) {

  # check color scale bin for each expression value
  bins <- c()
  for (i in 1:length(difexp)) {
    # which color bins are smaller than our difexp value
    # (for probet: i, mode:mode)
    inds <- which(difexp[[i]] > mybreaks)
    if (length(inds) == 0) {
      bins[[i]] <- 1
    } else if (length(inds) > 0)  {
      bins[[i]] <- max(inds) + 1
    }
  }	
  bins
}


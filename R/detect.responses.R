#
# Copyright (C) 2008-2011 Olli-Pekka Huovilainen and Leo Lahti 
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

  ######################################################################

  # Before I put a sketch on paper, the whole idea is worked out
  # mentally. In my mind I change the construction, make improvements,
  # and even operate the device. Without ever having drawn a sketch I
  # can give the measurements of all parts to workmen, and when
  # completed all these parts will fit, just as certainly as though I
  # had made the actual drawings. It is immaterial to me whether I run
  # my machine in my mind or test it in my shop. The inventions I have
  # conceived in this way have always worked. In thirty years there
  # has not been a single exception. My first electric motor, the
  # vacuum wireless light, my turbine engine and many other devices
  # have all been developed in exactly this way.
  #
  #                                                    - Nicola Tesla

  ######################################################################

detect.responses <-
function(datamatrix,
         network,
         initial.responses = 1,   # initial number of components. FIXME: is this used?
         max.responses = 10,      # max. responses # FIXME: check if really used
         max.subnet.size = 10,    # max. subnetwork size
         verbose = TRUE,          # print proc. information
         prior.alpha    = 1,      # Prior parameters
         prior.alphaKsi = 0.01,   # for VDP mixture
         prior.betaKsi  = 0.01,   # scale parameter for inverse Gamma
         update.hyperparams = 0,  # update hyperparameters. FIXME: check if this is applicable.
         implicit.noise = 0,      # Add implicit noise in vdp.mk.log.lambda.so and vdp.mk.hp.posterior.so 
         vdp.threshold = 1.0e-5,  # min. free energy improvement that stops VDP
         merging.threshold = 0,   # min. cost improvement for merging
         ite = Inf,                # max. iterations in updatePosterior
         information.criterion = "BIC", # information criterion for model selection
         speedup = TRUE,                 # speed up calculations by approximations
         speedup.max.edges = 10  # max. new joint models to be calculated; MI-based prefiltering applied
         )

{

#  NetResponse: Detect condition-specific network responses, given
#  network and a set of measurements of node activity in a set of
#  conditions.
#          
#  res <- netresponse(dataset, network)
#  
#  INPUT:
#
#  datamatrix: 'samples x features' matrix
#  network: binary matrix defining network between the features
#
#  OUTPUT:
#  Returns a set of subnetworks and their estimated
#  context-specific responses.
# 
#  res.costs are the cost function values at each state
#  res.moves has the indices of groups joined at each state in its columns
#  res.groupings holds the groupings at each level of the hierarchy
#  res.models has compressed representations of the models from each step

  set.seed(2341)
#  require(Matrix)

  
  # store here all params used in the model (defined in function call)
  params <- list(initial.responses = initial.responses, 
  	    	 max.responses = max.responses,
		 max.subnet.size = max.subnet.size,
                 verbose = verbose, 
		 prior.alpha = prior.alpha, 
		 prior.alphaKsi = prior.alphaKsi,
                 prior.betaKsi = prior.betaKsi, 
		 update.hyperparams = update.hyperparams, 
		 implicit.noise = implicit.noise,
                 vdp.threshold = vdp.threshold,
                 merging.threshold = merging.threshold,
		 ite = ite, information.criterion = information.criterion
		 )

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

  if (verbose) message("matching the features between network and datamatrix")  
  samples <- rownames(datamatrix)
  datamatrix   <- matrix(datamatrix[, network.nodes], nrow(datamatrix))
  colnames(datamatrix) <- network.nodes
  rownames(datamatrix) <- samples
  rm(samples)


  
#################################################################################

### INITIALIZE ###
  
  Nlog <- log( nrow( datamatrix ) )
  dim0 <- ncol( datamatrix )
  C    <- 0

  # print(" diagonal contains number of parameters for corresponding model")
  # off-diagonal tells number of parameters for the joint model
  # initial costs for the independent and joint models

  # Each variable pair has cost function value; add this as network edge property
  # is sparse matrices are slow, use arrays; they are faster than matrices
  #delta     <- matrix(Inf, dim, dim)
  #delta.nodes <- vector(Inf, network.nodes)

  # Network rows:                                        
  # 1) nodes 2) nodes 3) merging cost function delta
  rownames(network) <- c("node1", "node2")
  delta <- rep(NA, ncol(network))
  # FIXME: add individual models separately                                        
    
  # Storage list for calculated models
  model.nodes <- vector(length = length(network.nodes), mode = "list" ) # individual nodes
  model.pairs <- vector(length = ncol(network), mode = "list" ) # model for each pair

  gc()
  
########################################################################

### INDEPENDENT MODEL FOR EACH VARIABLE ###
  
  if (verbose) { cat("Compute cost for each variable\n") }

  for (k in 1:length(network.nodes)){

    node <- network.nodes[[k]]
    
    if ( verbose ) { cat(paste('Computing model for node', k, "/", ncol( datamatrix ), '\n')) }

    model <- vdp.mixt( matrix(datamatrix[, node], nrow( datamatrix )),
                      implicit.noise = implicit.noise,
                      prior.alpha = prior.alpha,
                      prior.alphaKsi = prior.alphaKsi,
                      prior.betaKsi = prior.betaKsi,
                      threshold = vdp.threshold,
                      initial.K = initial.responses,
                      ite = ite,
                      c.max = max.responses - 1,
                      speedup = speedup )

    cost.ind <- information.criterion(model$posterior$Nparams, Nlog, -model$free.energy, criterion = information.criterion) # COST for model
    C              <- C + cost.ind # Total cost
    model.nodes[[k]] <- pick.model.parameters(model, node)

}

gc()
  
if (verbose) { cat('done\n') }

  
##########################################################################################

###   compute costs for combined variable pairs  ###

for (edge in 1:ncol(network)){

  if (verbose) { cat(paste('Computing delta values for edge ', edge, '/', ncol(network), '\n')) }

  a <- network[1, edge]
  b <- network[2, edge]  
  vars            <- network.nodes[c(a, b)]
  model           <- vdp.mixt(
                              matrix(datamatrix[, vars], nrow( datamatrix )),
                              implicit.noise = implicit.noise,
                              prior.alpha = prior.alpha,
                              prior.alphaKsi = prior.alphaKsi,
                              prior.betaKsi = prior.betaKsi,
                              threshold = vdp.threshold,
                              initial.K = initial.responses,
                              ite = ite,
                              c.max = max.responses - 1,
                              speedup = speedup)

  # Compute COST-value for two independent subnets vs. joint model 
  # Negative free energy (-cost) is (variational) lower bound for P(D|H)
  # Use it as an approximation for P(D|H)
  # Cost for the indpendent and joint models
  # -cost is sum of two independent models (cost: appr. log-likelihoods)
  costind.ab     <-  information.criterion(model.nodes[[a]]$Nparams + model.nodes[[b]]$Nparams, Nlog, -(model.nodes[[a]]$free.energy + model.nodes[[b]]$free.energy), criterion = information.criterion)
  costjoint.ab   <-  information.criterion(model$posterior$Nparams, Nlog, -model$free.energy, criterion = information.criterion)
 
      # NOTE: COST is additive so summing is ok
      # change (increase) of the total COST / cost
  delta[[edge]] <- as.numeric(costjoint.ab - costind.ab)
      # Store these only if it would improve the cost; otherwise never needed
  if (-delta[[edge]] > merging.threshold) {    
    model.pairs[[edge]] <- pick.model.parameters(model, vars)
  } else {
    model.pairs[[edge]] <- 0
  }      
}

gc()

#######################################################################################


### MERGE VARIABLES ###

G <- lapply(1:ncol( datamatrix ), function( x ){ x }) # Place each node in a singleton subnet
move.cost.hist  <- matrix(c(0, 0, C), nrow = 3)

# if there are groups left sharing a link and improvement (there are
# connected items that have delta<0) then continue merging
# note that diag(network) has been set to 0
while ( !is.null(network) && any( -delta > merging.threshold )){

  if ( verbose ) { cat(paste('Combining groups, ', sum(!is.na(G)), ' group(s) left...\n'))} else{}
    
    # Identify the best neighbor pair in the network (also check that
    # the new merged pair would not exceed the max allowed subnetwork
    # size)

  tmp <- find.best.neighbor3(G, max.subnet.size, network, delta)
  delta <- tmp$delta

  # If merging still possible
  if (-tmp$mindelta > merging.threshold) {
    a <- tmp$a 
    b <- tmp$b
    best.edge <- tmp$best.edge
    # Store results                                        
    C <- C + tmp$mindelta
    move.cost.hist <- cbind(move.cost.hist, matrix(c(a, b, C), 3))

    # put the new group to a's place only for those variables for
    # which this is needed.  For others, put Inf on the a neighborgs,

    # combine a and b in the network, remove self-link a-a, remove b (row and col)
    tmp.join <- join.subnets2(network, delta, best.edge)
    network <- tmp.join$network
      delta <- tmp.join$delta
    model.nodes[[a]] <- model.pairs[[best.edge]]
    model.nodes[[b]] <- NA
    
    # remove self-links
    #keep <- (!network[1,] == network[2,])
    keep <- !(network[1,] == network[2,])
    network <- network[, keep]
    delta <- delta[keep]
    model.pairs <- model.pairs[keep]    

    # Merge groups G[[a]], G[[b]]
    G[[a]] <- sort(c(G[[a]], G[[b]]))
    G[[b]] <- NA
 
    # Skip the first b-1 elements as we only apply lower triangle here
    if ( ncol(network) == 1 ) {
      if ( verbose ) { cat("All nodes have been merged.\n") }
      delta <- Inf #indicating that no merging can be be done any more
    } else {
      # Compute new joint models for the new merged subnet and its neighborghs
      merge.edges <- which(is.na(delta))

      if (speedup && length(merge.edges) > speedup.max.edges) {
        # To speed up computation, pre-filter the edge set for which
        # new models are calculated.  Calculate empirical mutual
        # information between the first principal components of each
        # subnetwork pair. If number of new subnetwork pairs exceeds
        # the threshold, then calculate new model only for the
        # subnetwork pairs that have the highest mutual information.
        # It is expected that the subnetwork pair that will benefit
        # most from joint modeling will also be among the top mutual
        # infomation candidates. This way we can avoid calculating
        # exhaustive many models on large network hubs at each
        # update.
        #require(minet)
        mis <- c()
        mi.cnt <- 0  
        for (edge in which(is.na(delta))){
          mi.cnt <- mi.cnt + 1
          # Pick node indices
          a <- network[1, edge]
          i <- network[2, edge]
          dat <- cbind(prcomp(matrix(datamatrix[, network.nodes[G[[a]]]], nrow(datamatrix)), center = TRUE)$x,
                       prcomp(matrix(datamatrix[, network.nodes[G[[i]]]], nrow(datamatrix)), center = TRUE)$x)
          mis[[mi.cnt]] <- build.mim(dat, estimator="mi.empirical", disc = "equalwidth")[1,2]
        }
        merge.edges <- which(is.na(delta))[order(mis, decreasing = TRUE)[1:speedup.max.edges]]
        other.edges <- setdiff(which(is.na(delta)), merge.edges)
        delta[other.edges] <- Inf # needs Inf although not calculated; NAs would be confused with other merges later; models to be calculated are taken from is.na(delta) at each step so we cannot leave NAs there
      }
      
      for (edge in merge.edges){

        # Pick node indices
        a <- network[1, edge]
        i <- network[2, edge]
        vars  <- network.nodes[sort(c(G[[a]], G[[i]]))]

        model <- vdp.mixt(matrix(datamatrix[, vars], nrow( datamatrix )),
                          implicit.noise = 0,
                          prior.alpha = prior.alpha,
                          prior.alphaKsi = prior.alphaKsi,
                          prior.betaKsi = prior.betaKsi,
                          threshold = vdp.threshold,
                          initial.K = initial.responses,
                          ite = ite,
                          c.max = max.responses - 1,
                          speedup = speedup)

          # Store the joint models
          # cost for two independent vs. joint model
          # Negative free energy is (variational) lower bound for P(D|H)
          # Use this to approximate P(D|H)
        if (is.finite(model$free.energy)) {
          cost.ind <- information.criterion((model.nodes[[a]]$Nparams + model.nodes[[i]]$Nparams), Nlog, -(model.nodes[[a]]$free.energy + model.nodes[[i]]$free.energy), criterion = information.criterion)
          cost.joint <- information.criterion(model$posterior$Nparams, Nlog, -model$free.energy, criterion = information.criterion)
          
          # change (increase) of the total cost
          delta[[edge]] <- cost.joint - cost.ind
        } else  {
          warning("No free energy obtained.")
          delta[[edge]] <- Inf
        }
            
          if (-delta[[edge]] > merging.threshold) {  
            # Store joint model only if it would improve the cost
            model.pairs[[edge]] <- pick.model.parameters(model, vars)
          } else {
            model.pairs[[edge]] <- 0
          }
      }
    }
  } else{
    if ( verbose ) {cat(paste('Merging completed: no groups having links any more, or cost function improvement does not exceed the threshold.\n'))}
    break
  }
}
  
  # Remove left-out nodes (from the merges)
  nainds <- is.na(model.nodes)
  model.nodes <- model.nodes[!nainds]
  G <- G[!nainds]

  # Form a list of subnetworks (no filters)
  subnet.list <- lapply(G, function(x) { network.nodes[unlist(x)] })

  # name the subnetworks
  names(model.nodes) <- names(subnet.list) <- names(G) <- paste("Subnet-", 1:length(G), sep = "")  

  gc()

  # Convert original network to graphNEL (not before, to save more memory for computation stage)
  network.orig <- igraph.to.graphNEL(graph.data.frame(as.data.frame(t(network.orig)), directed = FALSE, vertices = data.frame(cbind(1:length(network.nodes), network.nodes))))
  nodes(network.orig) <- network.nodes
  
  new("NetResponseModel",
      moves = matrix(move.cost.hist, 3),
      last.grouping = G,     # network nodes given in indices
      subnets = subnet.list, # network nodes given in feature names; FIXME: remove available from models and G
      params = params,
      datamatrix = datamatrix,
      network = network.orig,
      models = model.nodes
      )
}

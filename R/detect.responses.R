#
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
         max.responses = 10,      
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
         information.criterion = "AIC", # information criterion for model selection
         speedup = TRUE,                 # speed up calculations by approximations
         speedup.max.edges = 10,  # max. new joint models to be calculated; MI-based prefiltering applied
	 mc.cores = 1, # number of cores for parallelization
         ... # Further arguments
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

  datamatrix <- check.matrix(datamatrix)

  tmp <- check.network(network, datamatrix, verbose = verbose)
  network <- tmp$formatted
  network.orig <- tmp$original
  delta <- tmp$delta
  network.nodes <- tmp$nodes
  rm(tmp)
  
  ### INITIALIZE ###
  
  if (verbose) message("matching the features between network and datamatrix")  
  samples <- rownames(datamatrix)
  datamatrix   <- matrix(datamatrix[, network.nodes], nrow(datamatrix))
  colnames(datamatrix) <- network.nodes
  rownames(datamatrix) <- samples
  rm(samples)

  Nlog  <- log( nrow( datamatrix ) ) # FIXME move this to params from all places
  nbins <- floor(sqrt(nrow(datamatrix))) # FIXME move this to params from all places

  # Store here all params used in the model (defined in function call)
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
		 ite = ite, 
		 information.criterion = information.criterion,
		 speedup = speedup,
		 speedup.max.edges = speedup.max.edges,
		 Nlog = Nlog,
		 nbins = nbins,
		 mc.cores = mc.cores
		 )

  # Place each node in a singleton subnet
  G <- lapply(1:ncol( datamatrix ), function( x ){ x }) 

  tmp <- filter.netw(network, delta, datamatrix, params)
  network <- tmp$network      
  delta <- tmp$delta

  gc()

  ########################################################################

  ### INDEPENDENT MODEL FOR EACH VARIABLE ###
  tmp <- independent.models(datamatrix, params)
  model.nodes <- tmp$nodes
  C <- sum(tmp$C)

  ###   compute costs for combined variable pairs  ###
  tmp <- pick.model.pairs(network, network.nodes, model.nodes, datamatrix, params) 
  model.pairs <- tmp$model.pairs
  delta <- tmp$delta

  #######################################################################################

  ### MERGE VARIABLES ###

  move.cost.hist  <- matrix(c(0, 0, C), nrow = 3)

  if (params$max.subnet.size > 1) {

    # if there are groups left sharing a link and improvement (there are
    # connected items that have delta<0) then continue merging
    # note that diag(network) has been set to 0
    while ( !is.null(network) && any( na.omit(-delta) > merging.threshold )){

      if ( verbose ) { message(paste('Combining groups, ', sum(!is.na(G)), ' group(s) left...\n'))} else{}
    
      # Identify the best neighbor pair in the network (also check that
      # the new merged pair would not exceed the max allowed subnetwork
      # size)

      tmp <- find.best.neighbor(G, max.subnet.size, network, delta)
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

        tmp.join <- join.subnets(network, delta, best.edge)
        network <- tmp.join$network
        delta <- tmp.join$delta
        model.nodes[[a]] <- model.pairs[[best.edge]]
        model.nodes[[b]] <- NA
    
        # remove self-links
        keep <- !(network[1,] == network[2,])
        network <- network[, keep]
        delta <- delta[keep]
        model.pairs <- model.pairs[keep]    

        # Merge groups G[[a]], G[[b]]
        G[[a]] <- sort(c(G[[a]], G[[b]]))
        G[[b]] <- NA

        # Skip the first b-1 elements as we only apply lower triangle here
        if ( ncol(network) == 1 ) {
          if ( verbose ) { message("All nodes have been merged.\n") }
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
            merge.edges <- which(is.na(delta))[order(get.mis(datamatrix, network, delta, network.nodes, G, params), decreasing = TRUE)[1:speedup.max.edges]]

	    # needs Inf although not calculated; NAs would be confused with other merges later; 
	    # models to be calculated are taken from is.na(delta) at each step so we cannot leave NAs there
            delta[setdiff(which(is.na(delta)), merge.edges)] <- Inf 

          }
      
          # FIXME: lapply, or parallelize to speed up
          for (edge in merge.edges) {
	    tmp <- update.model.pair(datamatrix, delta, network, edge, network.nodes, G, params, model.nodes, model.pairs)
	    model.pairs <- tmp$model.pairs
	    delta <- tmp$delta
	  }
       }

  } else {
    if ( verbose ) { message(paste('Merging completed: no groups having links any more, or cost function improvement does not exceed the threshold.')) }
    break
  }
 }
}

  
  # Remove left-out nodes (from the merges)
  nainds <- is.na(model.nodes)
  model.nodes <- model.nodes[!nainds]
  G <- G[!nainds]

  # Form a list of subnetworks (no filters)
  subnet.list <- lapply(G, function(x) { network.nodes[unlist(x)] }) # mclapply was slower here

  # name the subnetworks
  names(model.nodes) <- names(subnet.list) <- names(G) <- paste("Subnet-", 1:length(G), sep = "")  
  gc()

  # Convert original network to graphNEL (not before, to save more memory for computation stage)
  network.orig <- igraph.to.graphNEL(graph.data.frame(as.data.frame(t(network.orig)), directed = FALSE, vertices = data.frame(cbind(1:length(network.nodes), network.nodes))))
  nodes(network.orig) <- network.nodes

  
  model <- new("NetResponseModel",
      moves = matrix(move.cost.hist, 3),
      last.grouping = G,     # network nodes given in indices
      subnets = subnet.list, # network nodes given in feature names; FIXME: remove available from models and G
      params = params,
      datamatrix = datamatrix,
      network = network.orig,
      models = model.nodes
      )

}

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
         ite = Inf                # max. iterations in updatePosterior
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
		 ite = ite 
		 )

  # FIXME: later add other forms of sparse matrices from Matrix package
  accepted.formats <- c("matrix", "Matrix", "data.frame", "dgCMatrix", "dgeMatrix")
 
  # ensure datamatrix is a matrix
  if (!is.matrix(datamatrix)) {
    if (class(datamatrix) %in% accepted.formats) {
      message("Converting the input data into (sparse) matrix format.")
      datamatrix <- matrix(datamatrix)
    } else {
     stop(paste("datamatrix needs to be in one of the following formats:", paste(accepted.formats, collapse = "; ")))
    }    
  }

  # ensure network is a matrix
  if (!is.matrix(network)) {
    if (class(network)[[1]] %in% accepted.formats) {
#      message("Converting the input network into (sparse) matrix format.")
#      network <- Matrix(network)
#       network <- as.matrix(network) 
      # FIXME: dtpMatrix Triangular real matrices in packed storage (triangle only)
      # would save even more space than using general sparse real matrix here. Change.
      
    } else {
      stop(paste("network needs to be in one of the following formats:", paste(accepted.formats, collapse = "; ")))
    }
  }

  # remove self-links i.e. set network diagonal to zero
  network <- as.matrix(network)
  diag( network ) <- 0
  #network <- Matrix(network)

  # match the features between network and datamatrix
  # the names need to match if names are given
  if ( !is.null(rownames( network )) && !is.null(rownames( datamatrix )) ) {
    # pick samples that are in both network and response matrix                                        
    common.feats <- intersect(rownames(network), colnames(datamatrix))
    network      <- network[common.feats, common.feats]
    net.datamatrix   <- matrix(datamatrix[, common.feats], nrow(datamatrix))
    rownames(net.datamatrix) <- rownames(datamatrix)
    colnames(net.datamatrix) <- common.feats
    datamatrix <- net.datamatrix
    rm(net.datamatrix)
    
  } else if ( !nrow(network) == ncol(network) ) {

    stop("Error: symmetric network required.\n")

  } else if ( !nrow(network) == ncol(datamatrix) ) {

    stop("Error: Equal number of features required for the network and data matrix when feature names not given.\n")

  } else { warning("Warning: network and/or data features not named; matched by order.\n") }

  # store the original network (self-links removed, ordered to match the datamatrix)
  #  network <- Matrix(network)
  network.orig <- network

#################################################################################

### INITIALIZE ###
  
  Nlog <- log( nrow( datamatrix ) )
  dim0 <- dim <- ncol( datamatrix )
  C    <- 0

  # print(" diagonal contains number of parameters for corresponding model")
  # off-diagonal tells number of parameters for the joint model
  # initial costs for the independent and joint models

  # is sparse matrices are slow, use arrays; they are faster than matrices
  delta     <- matrix(Inf, dim, dim)

  # Storage list for calculated models
  model.list <- list()
  for (i in 1:nrow(network)) {  
    model.list[[i]] <- vector(length = i, "list")
  }
  # nested lists; the first level corresponds to rows and second level to cols of the other matrices
  # just lower-diagonal used

  gc()
  
########################################################################

### INDEPENDENT MODEL FOR EACH VARIABLE ###
  
  if (verbose) { cat("Compute cost for each variable\n") }

  for (k in 1:dim){

    if ( verbose ) {cat(paste('Computing model for variable', k, "/", dim, '\n'))}    

    model           <- vdp.mixt( matrix(datamatrix[, k], nrow( datamatrix )),
                                implicit.noise = implicit.noise,
                                prior.alpha = prior.alpha,
                                prior.alphaKsi = prior.alphaKsi,
                                prior.betaKsi = prior.betaKsi,
                                threshold = vdp.threshold,
                                initial.K = initial.responses,
                                ite = ite,
                                c.max = max.responses - 1 )
    # Save memory; include these later if needed
    model$prior <- model$opts <- NULL
    bic.ind <- bic(model$posterior$Nparams, Nlog, -model$free.energy) # BIC for model
    C              <- C + bic.ind # Total cost
    model.list[[k]][[k]] <- pick.model.parameters(model, colnames(datamatrix)[[k]])

}

gc()
  
if (verbose) { cat('done\n') }

##########################################################################################

###   compute costs for combined variable pairs  ###

for (a in 1:(dim - 1)){

  if (verbose) { cat(paste('Computing delta values for variable ', a, '/', dim, '\n')) }

  for (b in (a + 1):dim){

    # Require that the combined groups are connected in the network
    if (network[a, b]){ 

      vars            <- c(a, b)
      model           <- vdp.mixt(
                                  matrix(datamatrix[, vars], nrow( datamatrix )),
                                  implicit.noise = implicit.noise,
                                  prior.alpha = prior.alpha,
                                  prior.alphaKsi = prior.alphaKsi,
                                  prior.betaKsi = prior.betaKsi,
                                  threshold = vdp.threshold,
                                  initial.K = initial.responses,
                                  ite = ite,
                                  c.max = max.responses - 1)
      # Save memory; include these later if needed
      model$prior <- model$opts <- NULL

      # Compute BIC-value for two independent subnets vs. joint model 
      # Negative free energy (-cost) is (variational) lower bound for P(D|H)
      # Use it as an approximation for P(D|H)
      # Cost for the indpendent and joint models
      # -cost is sum of two independent models (cost: appr. log-likelihoods)
      bicind.ab     <-  bic(model.list[[a]][[a]]$Nparams + model.list[[b]][[b]]$Nparams, Nlog, -(model.list[[a]][[a]]$free.energy + model.list[[b]][[b]]$free.energy))
      bicjoint.ab   <-  bic(model$posterior$Nparams, Nlog, -model$free.energy)
      
      # NOTE: BIC is additive so summing is ok
      # change (increase) of the total BIC / cost
      dab <- bicjoint.ab - bicind.ab

      # Store these only if it would improve the cost; otherwise never needed
      if (-dab > merging.threshold) {
        #Nparams[a, b] <- model$posterior$Nparams # number of parameters in joint model
        delta[a, b] <- dab
        model.list[[a]][[b]] <- pick.model.parameters(model, colnames(datamatrix)[vars])
      } else {
        model.list[[a]][[b]] <- NA
      }      
    }
  }
}

gc()
  
#######################################################################################

### MERGE VARIABLES ###

G <- lapply(1:dim, function( x ){ x }) # Place each node in a singleton subnet
move.cost.hist  <- matrix(c(0, 0, C), nrow = 3)

for (j in 2:dim0){

  # if there are groups left sharing a link and improvement (there are
  # connected items that have delta<0) then continue merging
  # note that diag(network) has been set to 0
  if ( sum(network) > 0 && any( -delta > merging.threshold )){

    if ( verbose ) { cat(paste('Combining groups, ', nrow(network) - 1, ' group(s) left...\n'))} else{}
    gc() # clean up memory
    
    # Identify the best neighbor pair in the network (also check that
    # the new merged pair would not exceed the max allowed subnetwork
    # size)
    tmp <- find.best.neighbor(delta, G, max.subnet.size, network)
      a <- min(tmp$a, tmp$b)
      b <- max(tmp$a, tmp$b)

    # Store results
    C <- C + tmp$mindelta
    #   move.cost.hist <- cBind(move.cost.hist, Matrix(c(a, b, C), 3)) 
    move.cost.hist <- cbind(move.cost.hist, matrix(c(a, b, C), 3))     

    # put the new group to a's place only for those variables for
    # which this is needed.  For others, put Inf on the a neighborgs,
    
    # combine a and b in the network, remove self-link a-a, remove b (row and col)
    network <- join.subnets(network, a, b)
    model.list[[a]][[a]] <- model.list[[a]][[b]]

    # Merge groups G[[a]], G[[b]]
    G[[a]] <- sort(c(G[[a]], G[[b]]))
    G      <- G[-b]
    
    # remove the merged group
    # network b already removed in join.subnets() above
    delta     <- delta[-b, -b]

    # remove bth elements in model list
    model.list[[b]] <- NULL
    # if there are more rows available 
    # then remove the b:th element from each..
    for (b.idx in 1:length(model.list)) {    
      model.list[[b.idx]][[b]] <- NULL
    } 
    
    # Skip the first b-1 elements as we only apply lower triangle here
    if ( nrow(network) == 1 ) {
      if ( verbose ) { cat("All nodes have been merged.\n") }
      delta <- Inf #indicating that no merging can be be done any more
    } else {

      # Infinite joint costs etc with a for groups not linked to a
      # Note that for Nparams we need also a-a information    
      # Therefore do not replace a, a
      delta[a, ]     <- delta[, a]     <- Inf

      # Compute new joint models for a and its neighborghs
      for (i in 1:nrow(network)){

        # compute combined model only if a and i are linked
        if (network[a, i] & length(c(G[[a]], G[[i]])) <= max.subnet.size){
          vars  <- sort(c(G[[a]], G[[i]]))

          model <- vdp.mixt(matrix(datamatrix[, vars], nrow( datamatrix )),
                                implicit.noise = implicit.noise,
                                prior.alpha = prior.alpha,
                                prior.alphaKsi = prior.alphaKsi,
                                prior.betaKsi = prior.betaKsi,
                                threshold = vdp.threshold,
                                initial.K = initial.responses,
                                ite = ite,
                                c.max = max.responses - 1 )
          # Save memory; include these later if needed
          model$prior <- model$opts <- NULL

          # Store the joint models (always a < i)
          minind <- min(c(a,i))
	  maxind <- max(c(a,i))

          # BIC-cost for two independent vs. joint model
          # Negative free energy (-cost) is (variational) lower bound for P(D|H)
          # Use this to approximate P(D|H)

          big.ind   <- bic((model.list[[a]][[a]]$Nparams + model.list[[i]][[i]]$Nparams), Nlog, -(model.list[[a]][[a]]$free.energy + model.list[[i]][[i]]$free.energy))
          big.joint <- bic(model$posterior$Nparams, Nlog, -model$free.energy)
          
          # change (increase) of the total BIC (cost)
          delta[minind, maxind] <- big.joint - big.ind

          if (-delta[minind, maxind] > merging.threshold) {  
            # Store joint model only if it would improve the cost
            model.list[[minind]][[maxind]] <- pick.model.parameters(model, colnames(datamatrix)[vars])
          } else {
            model.list[[minind]][[maxind]] <- NA
          }
        }
      }
    }
  } else{
    if ( verbose ) {cat(paste('Merging completed: no groups having links any more, or no improvement possible on level', j, '\n'))}
    break
  }
}

  gc()

  #if (is.null(colnames(datamatrix))) { nodes <- as.character(1:ncol(datamatrix)) } else {nodes <- colnames(datamatrix)}
  #if (is.null(rownames(datamatrix))) { samples <- as.character(1:nrow(datamatrix)) } else {samples <- rownames(datamatrix)}

  if (is.null(colnames(datamatrix))) { colnames(datamatrix) <- as.character(1:ncol(datamatrix)) }
  if (is.null(rownames(datamatrix))) { rownames(datamatrix) <- as.character(1:nrow(datamatrix)) }  
    
  # Form a list of subnetworks (no filters)
  #subnet.list <- lapply(G, function(x) { nodes[unlist(x)] })
  subnet.list <- lapply(G, function(x) { rownames(network.orig)[unlist(x)] })  

  # Pick the final subnetwork models from the model.list (the diagonal)
  diag.models <- list()
  for (k in rev(1:length(model.list))) { # go through in reverse order to remove last k at each iter.
    diag.models[[k]] <- model.list[[k]][[k]]
    model.list[[k]] <- NULL # saving memory
  }
    
  # name the subnetworks
  names(diag.models) <- names(subnet.list) <- names(G) <- paste("Subnet-", 1:length(G), sep = "")  

  gc()
  
  new("NetResponseModel",
      moves = matrix(move.cost.hist, 3),
      last.grouping = G,     # network nodes given in indices
      subnets = subnet.list, # network nodes given in feature names; FIXME: remove available from models and G
      params = params,
      datamatrix = datamatrix,
      network = network.orig,
      models = diag.models
      )      

}

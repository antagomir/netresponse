detect.responses <-
function(datamatrix,
         network,
         initial.responses = 1,  # initial number of components. FIXME: is this used?
         max.responses = 10,     # max. responses # FIXME: check if really used
         max.subnet.size = 20,   # max. subnetwork size
         verbose = TRUE,         # print proc. information
         prior.alpha    = 1,     # Prior parameters
         prior.alphaKsi = 0.01,  # for VDP mixture
         prior.betaKsi  = 0.01,  # scale parameter for inverse Gamma
         update.hyperparams = 0, # update hyperparameters. FIXME: check if this is applicable.
         implicit.noise = 0, # Add implicit noise in vdp.mk.log.lambda.so and vdp.mk.hp.posterior.so 
         threshold = 1.0e-5, # min. free energy improvement that stops VDP
         ite = Inf          # max. iterations in updatePosterior
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
#
# Copyright (C) 2008-2010 Olli-Pekka Huovilainen and Leo Lahti
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

  set.seed(2341)
  
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
                 threshold = threshold, 
		 ite = ite 
		 )

  # FIXME: later add other forms of sparse matrices from Matrix package
  accepted.formats <- c("matrix", "Matrix", "data.frame", "dgCMatrix")
 
  # ensure datamatrix is a matrix
  if (!is.matrix(datamatrix)) {
    if (class(datamatrix) %in% accepted.formats) {
      message("Converting the input data into matrix format.")
      datamatrix <- as.matrix(datamatrix)
    } else {
     stop(paste("datamatrix needs to be in one of the following formats:", paste(accepted.formats, collapse = "; ")))
    }    
  }

  # ensure network is a matrix
  if (!is.matrix(network)) {
    if (class(network)[[1]] %in% accepted.formats) {
      message("Converting the input network into matrix format.")
      network <- as.matrix(network)
    } else {
      stop(paste("network needs to be in one of the following formats:", paste(accepted.formats, collapse = "; ")))
    }
  }
	      
  # network diagonal has to be zero, i.e. self-links not taken into account
  diag( network ) <- 0

  # match the features between network and datamatrix
  # the names need to match if names are given
  if ( !is.null(rownames( network )) && !is.null(rownames( datamatrix )) ) {
    # pick samples that are in both network and response matrix                                        
    common.feats <- intersect(rownames(network), colnames(datamatrix))
    network      <- network[common.feats, common.feats]
    net.datamatrix   <- matrix(datamatrix[, common.feats], nrow = nrow(datamatrix))
    rownames(net.datamatrix) <- rownames(datamatrix)
    colnames(net.datamatrix) <- common.feats
    datamatrix <- net.datamatrix

  } else if ( !nrow(network) == ncol(network) ) {

    stop("Error: symmetric network required.\n")

  } else if ( !nrow(network) == ncol(datamatrix) ) {

    stop("Error: Equal number of features required for the network and data matrix when feature names not given.\n")

  } else { warning("Warning: network and/or data features not named; matched by order.\n") }

  # store the original network (self-links removed, ordered to match the datamatrix)
  network.orig <- network


#################################################################################

### INITIALIZE ###
  
  Nlog <- log( nrow( datamatrix ) )
  dim0 <- dim <- ncol( datamatrix )
  C    <- 0

  # print(" diagonal contains number of parameters for corresponding model")
  # off-diagonal tells number of parameters for the joint model
  # initial costs for the independent and joint models
  #Nresponses <-  rep(NA, dim0)
  H         <-   rep(Inf, dim0) #array(Inf, dim = c(3, 1))
  costs     <- array(Inf, dim = c(dim, dim))
  Nparams   <- array(Inf, dim = c(dim, dim))
  bic.ind   <- array(Inf, dim = c(dim, dim))
  bic.joint <- array(Inf, dim = c(dim, dim))
  delta     <- array(Inf, dim = c(dim, dim))

  # Storage list for calculated models
  model.list <- list()
  for (i in 1:nrow(network)) {  
    model.list[[i]] <- vector(length = 1:i, "list")
  }
  # nested lists; the first level corresponds to rows and second level to cols of the other matrices
  # just lower-diagonal used

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
                                threshold = threshold,
                                initial.K = initial.responses,
                                ite = ite,
                                c.max = max.responses - 1 )
    
    H[[k]]         <- cost  <- model$free.energy # -cost is lower bound for log(P(D|H))
    Nparams[k, k]  <- model$posterior$Nparams # number of parameters for model
    bic            <- Nparams[k, k]*Nlog + 2*H[[k]] # BIC for model
    C              <- C + bic # Total cost, previously: C + H(k)

    model.list[[k]][[k]] <- pick.model.parameters(model, colnames(datamatrix)[[k]])

}

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
                                  threshold = threshold,
                                  initial.K = initial.responses,
                                  ite = ite,
                                  c.max = max.responses - 1)

      # Store the joint models
      model.list[[a]][[b]] <- pick.model.parameters(model, colnames(datamatrix)[vars])
    
      costs[a, b]     <- costs[b, a]   <- model$free.energy           # Store cost for joint model. 
      Nparams[a, b]   <- Nparams[b, a] <- model$posterior$Nparams     # number of parameters in joint model

      # Compute BIC-value for two independent subnets vs. joint model 
      # Negative free energy (-cost) is (variational) lower bound for P(D|H)
      # Use it as an approximation for P(D|H)
      # Cost for the indpendent and joint models
      # -cost is sum of two independent models (H: appr. log-likelihoods)
      bic.ind[a, b] <- bic.ind[b, a] <- (Nparams[a, a] + Nparams[b, b])*Nlog + 2*(H[[a]] + H[[b]])
      bic.joint[a, b] <- bic.joint[b, a] <- Nparams[a, b]*Nlog + 2*(costs[a, b]) 
      # = Nparams[a, b]*Nlog - 2*(-costs[a, b])

      # NOTE: BIC is additive so summing is ok
      # change (increase) of the total BIC / cost
      delta[a, b] <- delta[b, a] <- bic.joint[a, b] - bic.ind[a, b]     
    }
  }
}


#######################################################################################

### MERGE VARIABLES ###
  
  
G <- lapply(1:dim, function( x ){ x }) # Place each node in a singleton subnet
move.cost.hist  <- matrix(c(0, 0, C), nrow = 3)

for (j in 2:dim0){

  # if there are groups left sharing a link and improvement (there are
  # connected items that have delta<0) then continue merging
  # note that diag(network) has been set to 0
  if (sum(network) > 0 && any( delta < 0 )){

    if ( verbose ) { cat(paste('Combining groups, ', nrow(network) - 1, ' group(s) left...\n'))} else{}

    # Identify the best neighbor pair in the network (also check that
    # the new merged pair would not exceed the max allowed subnetwork
    # size)
    tmp <- find.best.neighbor(delta, G, max.subnet.size, network)
      a <- tmp$a
      b <- tmp$b

    # Store results
    C <- C + tmp$mindelta
    move.cost.hist <- cbind(move.cost.hist, matrix(c(a, b, C), nrow = 3)) 

    # put the new group to a's place only for those variables for
    # which this is needed.  For others, put Inf on the a neighborgs,
    H[[a]] <- costs[a, b]

    # combine a and b in the network, remove self-link a-a, remove b (row and col)
    network <- join.subnets(network, a, b)

    # number of parameters for the subnets
    Nparams[a, a] <- Nparams[a, b]
    bic.ind[a, a] <- bic.joint[a, b]

    model.list[[a]][[a]] <- model.list[[a]][[b]]

    # Merge groups G[[a]], G[[b]]
    removed.group.vars <- G[[b]]
    G      <- G[-b]
    G[[a]] <- sort(c(as.numeric(G[[a]]), removed.group.vars))

    # remove the merged group
    # network b already removed in join.subnets() above
    #Nresponses <- Nresponses[-b]
    H          <- H[-b]
    bic.ind    <- bic.ind[-b, -b]
    bic.joint  <- bic.joint[-b, -b]
    costs      <- costs[-b, -b]
    Nparams    <- Nparams[-b, -b]
    delta      <- delta[-b, -b]

    # remove bth elements in model list
    model.list[[b]] <- NULL
    #if (length(model.list) >= 1) {
      # if there are more rows available 
      # then remove the b:th element from each..
    for (b.idx in 1:length(model.list)) {    
      model.list[[b.idx]][[b]] <- NULL
    } 
    #}
    
    # Skip the first b-1 elements as we only apply lower triangle here

    if ( nrow(network) == 1 ) {
      if ( verbose ) {cat("All nodes have been merged.\n")}
      delta > Inf #indicating that no merging can be be done any more
    } else {
    
      # Infinite joint costs etc with a for groups not linked to a
      # Note that for Nparams we need also a-a information    
      Nparams[a, -a] <- Nparams[-a, a] <- Inf
      bic.ind[a, -a] <- bic.ind[-a, a] <- Inf
      costs[a, ]     <- costs[, a]     <- Inf
      bic.joint[a, ] <- bic.joint[, a] <- Inf
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
                                threshold = threshold,
                                initial.K = initial.responses,
                                ite = ite,
                                c.max = max.responses - 1 )

          # Store the joint models (always a < i)
          minind <- min(c(a,i))
	  maxind <- max(c(a,i))
	  #print(c(minind,maxind, length(model.list)))
	  
          model.list[[minind]][[maxind]] <- pick.model.parameters(model, colnames(datamatrix)[vars])

          costs[a,i]   <- costs[i, a]   <- model$free.energy  # cost for joint model
          Nparams[a,i] <- Nparams[i, a] <- model$posterior$Nparams  # number of parameters in joint model

          # BIC-cost for two independent vs. joint model
          # Negative free energy (-cost) is (variational) lower bound for P(D|H)
          # Use this to approximate P(D|H)
          bic.ind[a, i]   <- bic.ind[i, a]   <- (Nparams[a, a] + Nparams[i, i])*Nlog + 2*(H[a] + H[i])
          bic.joint[a, i] <- bic.joint[i, a] <- Nparams[a, i]*Nlog + 2*(costs[a, i])

          # change (increase) of the total BIC (cost)
          delta[a, i] <- delta[i, a] <- bic.joint[a, i] - bic.ind[a, i]
        }
      }
    }
  } else{
    if ( verbose ) {cat(paste('Merging completed: no groups having links any more, or no improvement possible on level', j, '\n'))}
    break
  }
}

  if ( length(bic.ind) > 1 ) {
    costs <- diag(bic.ind)
  } else {
    costs <- bic.ind
  }
  
  if (is.null(colnames(datamatrix))) { nodes <- as.character(1:ncol(datamatrix)) } else {nodes <- colnames(datamatrix)}
  if (is.null(rownames(datamatrix))) { samples <- as.character(1:nrow(datamatrix)) } else {samples <- rownames(datamatrix)}
    
  # Form a list of subnetworks (no filters)
  subnet.list <- lapply(G, function(x) { nodes[unlist(x)] })

  # Pick the final subnetwork models from the model.list (the diagonal)
  diag.models <- list()
  for (k in 1:length(model.list)) {
    diag.models[[k]] <- model.list[[k]][[k]]
  }
  model.list <- diag.models

  # name the subnetworks
  names(model.list) <- names(subnet.list) <- names(G) <- names(costs) <- paste("Subnet-", 1:length(G), sep = "")  
        
  model <- 
  new("NetResponseModel",
      moves = matrix(move.cost.hist[1:2,], 2),
      costs = costs,
      last.grouping = G, # network nodes given in indices
      subnets = subnet.list, # network nodes given in feature names
      params = params,
      nodes = nodes,
      samples = samples,
      datamatrix = datamatrix,
      network = network.orig,
      models = model.list
      )      
   
}

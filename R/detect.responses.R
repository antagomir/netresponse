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

#' detect.responses
#' 
#' Main function of the NetResponse algorithm. 
#' Detect condition-specific network responses, given
#' network and a set of measurements of node activity in a set of
#' conditions. Returns a set of subnetworks and their estimated
#' context-specific responses.
#'
#' @param datamatrix Matrix of samples x features. For example, gene expression
#'   matrix with conditions on the rows, and genes on the columns. The matrix
#'   contains same features than the 'network' object, characterizing the network
#'   states across the different samples.
#' @param network Binary network describing undirected pairwise interactions between
#'   features of 'datamatrix'. The following formats are supported: binary
#'   matrix, graphNEL, igraph, graphAM, Matrix, dgCMatrix, dgeMatrix
#' @param initial.responses Initial number of components for each subnetwork
#'   model. Used to initialize calculations.
#' @param max.responses Maximum number of responses for each subnetwork. Can be
#'   used to limit the potential number of network states.
#' @param max.subnet.size Numeric. Maximum allowed subnetwork size.
#' @param verbose Logical. Verbose parameter.
#' @param implicit.noise Implicit noise parameter. Add implicit noise to vdp
#'   mixture model. Can help to avoid overfitting to local optima, if this
#'   appears to be a problem.
#' @param update.hyperparams Logical. Indicate whether to update
#'   hyperparameters during modeling.
#' @param prior.alpha,prior.alphaKsi,prior.betaKsi Prior parameters for
#'   Gaussian mixture model that is calculated for each subnetwork
#'   (normal-inverse-Gamma prior). alpha tunes the mean; alphaKsi and betaKsi are
#'   the shape and scale parameters of the inverse Gamma function, respectively.
#' @param vdp.threshold Minimal free energy improvement after which the
#'   variational Gaussian mixture algorithm is deemed converged.
#' @param merging.threshold Minimal cost value improvement required for merging
#' two subnetworks.
#' @param ite Defines maximum number of iterations on posterior update
#' (updatePosterior). Increasing this can potentially lead to more accurate
#' results, but computation may take longer.
#' @param information.criterion Information criterion for model selection.
#' Default is BIC (Bayesian Information Criterion); other options include AIC
#' and AICc.
#' @param speedup Takes advantage of approximations to PCA, mutual information
#' etc in various places to speed up calculations. Particularly useful with
#' large and densely connected networks and/or large sample size.
#' @param speedup.max.edges Used if speedup = TRUE. Applies prefiltering of
#' edges for calculating new joint models between subnetwork pairs when
#' potential cost changes (delta) are updated for a newly merged subnetwork and
#' its neighborghs. Empirical mutual information between each such subnetwork
#' pair is calculated based on their first principal components, and joint
#' models will be calculated only for the top candidates up to the number
#' specified by speedup.max.edges. It is expected that the subnetwork pair that
#' will benefit most from joint modeling will be among the top mutual
#' infomation candidates. This way it is possible to avoid calculating
#' exhaustive many models on the network hubs.
#' @param positive.edges Consider only the edges with positive association. Currently measured with Spearman correlation.
#' @param mc.cores Number of cores to be used in parallelization. See
#' help(mclapply) for details.
#' @param mixture.method Specify the approach to use in mixture modeling.
#' Options. vdp (nonparametric Variational Dirichlet process mixture model);
#' bic (based on Gaussian mixture modeling with EM, using BIC to select the
#' optimal number of components)
#' @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture with mixture.method = "bic"
#' @param pca.basis Transform data first onto PCA basis to try to avoid problems with non-diagonal covariances.
#' @param ... Further optional arguments to be passed.
#' @return NetResponseModel object.
#' @author Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse").
#' @keywords methods iteration
#' @export
#' @examples
#' data(toydata)        # Load toy data set
#' D    <- toydata$emat   # Response matrix (for example, gene expression)
#' netw <- toydata$netw   # Network
#' 
#' # Run NetReponse algorithm
#' model <- detect.responses(D, netw, verbose = FALSE)
detect.responses <- function(datamatrix,
         network = NULL,
         initial.responses = 1,   # initial number of components. FIXME: is this used?
         max.responses = 10,      
         max.subnet.size = 10,    # max. subnetwork size
         verbose = TRUE,          # print intermediate messages
         prior.alpha    = 1,      # Prior parameters
         prior.alphaKsi = 0.01,   # for VDP mixture
         prior.betaKsi  = 0.01,   # scale parameter for inverse Gamma
         update.hyperparams = 0,  # update hyperparameters. FIXME: check if this is applicable.
         implicit.noise = 0,      # Add implicit noise in vdp.mk.log.lambda.so and vdp.mk.hp.posterior.so 
         vdp.threshold = 1.0e-5,  # min. free energy improvement that stops VDP
         merging.threshold = 0,   # min. cost improvement for merging
         ite = Inf,               # max. iterations in updatePosterior
         information.criterion = "BIC", # information criterion for node merging
         speedup = TRUE,          # speed up calculations by approximations
         speedup.max.edges = 10,  # max. new joint models to be calculated; MI-based prefiltering applied
	 positive.edges = FALSE, # If TRUE, consider positive edges only	 
	 mc.cores = 1, # number of cores for parallelization
         mixture.method = "vdp", # Which approach to use for sample mixture estimation within given subnet. Options: bic/vdp
	 bic.threshold = 0,
	 pca.basis = FALSE,
	 ... # Further arguments
)

{

#fs <- list.files("~/Rpackages/netresponse/netresponse/R/", full.names = TRUE); for (f in fs) {source(f)}; datamatrix <- D; network <- netw; initial.responses = 1; max.responses = 3; max.subnet.size = 10; verbose = TRUE; prior.alpha = 1; prior.alphaKsi = 0.01; prior.betaKsi  = 0.01;	update.hyperparams = 0; implicit.noise = 0; vdp.threshold = 1.0e-5; merging.threshold = 1; ite = Inf; information.criterion = "BIC"; speedup = TRUE; speedup.max.edges = 10; mc.cores = 1; mixture.method = "bic"; bic.threshold = 0; pca.basis = FALSE          

  # Check data matrix validity         
  datamatrix <- check.matrix(datamatrix)

  # Check network validity and polish
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
		 Nlog = log( nrow( datamatrix ) ),
		 nbins = floor(sqrt(nrow(datamatrix))),
		 mc.cores = mc.cores,
		 mixture.method = mixture.method,
		 bic.threshold = bic.threshold,
		 positive.edges = positive.edges, 
		 pca.basis = pca.basis		 
		 )

  # Place each node in a singleton subnet
  G <- lapply(1:ncol( datamatrix ), function( x ){ x }) 

  # Filter network
  tmp <- filter.netw(network, delta, datamatrix, params)
  network <- tmp$network      
  delta <- tmp$delta
  # FIXME: for more efficient memory usage, remove from the datamatrix those nodes which are
  # not in the network. But check that the indices are not confused.

  gc()

  ########################################################################

  ### INDEPENDENT MODEL FOR EACH VARIABLE ###
  tmp <- independent.models(datamatrix, params)
  node.models <- tmp$nodes # model parameters
  C <- sum(tmp$C)

  ### MERGE VARIABLES ###

  # Store agglomeration steps
  move.cost.hist  <- matrix(c(0, 0, C), nrow = 3)

  if (params$max.subnet.size > 1) {

    ### compute costs for combined (singleton) variable pairs  ###
    tmp <- pick.model.pairs(network, network.nodes, node.models, datamatrix, params)

    model.pairs <- tmp$model.pairs
    delta <- tmp$delta

    # if there are groups left sharing a link and improvement (there are
    # connected items that have delta<0) then continue merging
    # note that diag(network) has been set to 0
    while ( !is.null(network) && any( na.omit(-delta) > merging.threshold )){

      if ( verbose ) { message(paste('Combining groups, ', sum(!is.na(G)), ' group(s) left...\n'))} else{}
    
      # Identify the best neighbor pair in the network (also check that
      # the new merged pair would not exceed the max allowed subnetwork
      # size)

      tmp <- find.best.neighbor(G, params$max.subnet.size, network, delta)
      delta <- tmp$delta
      best.edge <- tmp$best.edge

      # If merging is still possible
      if (-tmp$mindelta > merging.threshold) {

        a <- tmp$a 
        b <- tmp$b
        C <- C + tmp$mindelta
        move.cost.hist <- cbind(move.cost.hist, matrix(c(a, b, C), 3))      

        # put the new group to a's place only for those variables for
        # which this is needed.  For others, put Inf on the a neighborgs,
        # combine a and b in the network, remove self-link a-a, 
	# remove b (row and col)
        tmp.join <- join.subnets(network, delta, best.edge)
        network <- tmp.join$network
        delta <- tmp.join$delta
        node.models[[a]] <- model.pairs[[best.edge]]
        node.models[[b]] <- NA
    
        # remove self-links
        keep <- !(network[1,] == network[2,])
        network <- network[, keep]
        delta <- delta[keep]
        model.pairs <- model.pairs[keep]    

        # Merge groups G[[a]], G[[b]]
        G[[a]] <- sort(c(G[[a]], G[[b]]))
        G[[b]] <- NA

        # Skip the first b-1 elements as we only apply lower triangle here
        if ( ncol(network) <= 1 ) {
          if ( verbose ) { message("All nodes have been merged.\n") }
          delta <- Inf # indicates no merging can be be done any more
        } else {
          # Compute new joint models for the new merged subnet and its neighborghs
          merge.edges <- which(is.na(delta))
	   
	  # Remove edges that would exceed max.size
	  # FIXME: include as part of cost function?

	  if (length(merge.edges) > 0) {

  	    new.sizes <- apply(matrix(network[, merge.edges], 2), 2, function (x) {length(c(G[[x[[1]]]], G[[x[[2]]]]))})
	    merge.edges <- merge.edges[new.sizes <= params$max.subnet.size]

            if (speedup && length(merge.edges) > speedup.max.edges) {

              # To speed up computation, pre-filter the edge set for which
              # new models are calculated.  Calculate empirical associations
              # between the first principal components of each
              # subnetwork pair. If number of new subnetwork pairs exceeds
              # the threshold, then calculate new model only for the
              # subnetwork pairs that have the highest associations.
              # It is expected that the subnetwork pair that will benefit
              # most from joint modeling will also be among the top 
              # candidates. This way we can avoid calculating
              # exhaustive many models on large network hubs at each
              # update.
              merge.edges <- which(is.na(delta))[order(get.mis(datamatrix, network, delta, network.nodes, G, params), decreasing = TRUE)]
	    
	      # Remove edges that would exceed max.size
	      new.sizes <- apply(matrix(network[, merge.edges], 2), 2, function (x) {length(c(G[[x[[1]]]], G[[x[[2]]]]))})
	      keep <- (new.sizes <= params$max.subnet.size)
	      merge.edges <- merge.edges[keep][1:speedup.max.edges]
	    
	      # Needs Inf: NAs would be confused with other merges later since
	      # models to be calculated are taken from is.na(delta) at each step
              delta[setdiff(which(is.na(delta)), merge.edges)] <- Inf 

            } 

          }
      
          # TODO: parallelize to speed up
          for (edge in merge.edges) {
	    tmp <- update.model.pair(datamatrix, delta, network, edge, network.nodes, G, params, node.models, model.pairs)
	    model.pairs <- tmp$model.pairs
	    delta <- tmp$delta
	  }

	  # Didn't work straight away
	  #tmp <- mclapply(merge.edges, function (edge) {
	  #  update.model.pair(datamatrix, delta, network, edge, network.nodes, G, params, node.models, model.pairs)
	  #}, mc.cores = params$mc.cores)	  
	  #model.pairs <- lapply(tmp, function (x) {x$model.pairs})
	  #delta <- sapply(tmp, function (x) {x$delta})

       }

      } else {
        if ( verbose ) { message(paste('Merging completed: no groups having links any more, or cost function improvement does not exceed the threshold.')) }
        break
      }
    } 
  }
  
  # Remove left-out nodes (from the merges)
  nainds <- is.na(node.models)
  node.models <- node.models[!nainds]
  G <- G[!nainds]

  # Form a list of subnetworks (no filters)
  # mclapply was slower here  
  subnet.list <- lapply(G, function(x) { network.nodes[unlist(x)] }) 

  # name the subnetworks
  names(node.models) <- names(subnet.list) <- names(G) <- paste("Subnet-", 1:length(G), sep = "")  
  gc()

  # Convert original network to graphNEL (not before, to save more memory for computation stage)
  network.orig <- igraph.to.graphNEL(graph.data.frame(as.data.frame(t(network.orig)), directed = FALSE, vertices = data.frame(cbind(1:length(network.nodes), network.nodes))))
  nodes(network.orig) <- network.nodes

  # For one-dimensional subnets, 
  # order the modes by magnitude to simplify interpretation
  for (mi in 1:length(node.models)) {
    if (ncol(node.models[[mi]]$mu) == 1 && length(node.models[[mi]]$w) > 1) {
      o <- order(node.models[[mi]]$mu[,1])
      node.models[[mi]]$mu <- matrix(node.models[[mi]]$mu[o,], nrow = length(o))
      node.models[[mi]]$sd <- matrix(node.models[[mi]]$sd[o,], nrow = length(o))
      node.models[[mi]]$w <- node.models[[mi]]$w[o]
      rownames(node.models[[mi]]$mu) <- rownames(node.models[[mi]]$sd) <- names(node.models[[mi]]$w) <- paste("Mode-", 1:length(node.models[[mi]]$w), sep = "")
    }
  }

  # FIXME: if all nodes will be combined (merging.threshold = -Inf), there will be an error. Fix.
  #  costs: cost function values at each state
  #  moves: indices of groups joined at each state in its columns
  #  groupings: groupings at each level of the hierarchy
  #  models: compressed representations of the models from each step

  model <- new("NetResponseModel",
      moves = matrix(move.cost.hist, 3),
      last.grouping = G,     # network nodes given in indices
      subnets = subnet.list, # network nodes given in feature names; FIXME: remove available from models and G
      params = params,
      datamatrix = datamatrix,
      network = network.orig,
      models = node.models
      )
}


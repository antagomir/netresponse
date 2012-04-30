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
#' @return 
#'   \item{nodes }{Model for each node} 
#'   \item{C }{Costs for individual models}
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @export
#' @examples #
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
#' @examples #
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



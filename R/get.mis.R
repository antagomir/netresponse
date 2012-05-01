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
#' @author Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @examples #
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


edge.delta <- function (edge, network, network.nodes, datamatrix, params, model.nodes) {

    if ( params$verbose ) { message(paste('Computing delta values for edge ', edge, '/', ncol(network), '\n')) }
    a <- network[1, edge]
    b <- network[2, edge]
    vars <- network.nodes[c(a, b)]

    tmp <- mixture.model(matrix(datamatrix[, vars], nrow( datamatrix )), params$mixture.method, params$max.responses, params$implicit.noise, params$prior.alpha, params$prior.alphaKsi, params$prior.betaKsi, params$vdp.threshold, params$initial.responses, params$ite, params$speedup, params$bic.threshold) 
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
			
    return(list(model = mod.pair, delt = delt))
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



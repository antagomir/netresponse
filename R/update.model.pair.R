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
#' @examples #
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


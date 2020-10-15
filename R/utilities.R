#' Write NetResponse results summary into a file.
#' 
#' Experimental version.
#' 
#' @usage write.netresponse.results(x, subnet.ids = NULL, filename)
#' @param x NetResponseModel
#' @param subnet.ids List of subnet ids to consider. By default, all subnets.
#' @param filename Output file name.
#' @return Used for side effects.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation('netresponse')
#' @keywords utilities
write.netresponse.results <- function(x, subnet.ids = NULL, filename) {
    
    write("NetResponse - subnetworks", file = filename, append = FALSE)
    write("==========================\n", file = filename, append = TRUE)
    
    if (is.null(subnet.ids)) {
        subnet.ids <- names(x@subnets)
    }
    
    for (nam in subnet.ids) {
        write(nam, file = filename, append = TRUE)
        write(x@subnets[[nam]], file = filename, append = TRUE)
        write("----------------", file = filename, append = TRUE)
    }
    
}





#' get.model.parameters
#' 
#' Retrieve the mixture model parameters of the NetResponse algorithm for a
#' given subnetwork.
#' 
#' Only the non-empty components are returned. Note: the original data matrix
#' needs to be provided for function call separately.
#' 
#' @usage get.model.parameters(model, subnet.id = NULL)
#' @param model Result from NetResponse (detect.responses function).
#' @param subnet.id Subnet identifier. A natural number which specifies one of
#' the subnetworks within the 'model' object.
#' @return A list with the following elements: \item{mu}{ Centroids for the
#' mixture components. Components x nodes.} \item{sd}{ Standard deviations for
#' the mixture components. A vector over the nodes for each component, implying
#' the diagonal covariance matrix of the model (i.e. diag(std^2)). Components x
#' nodes} \item{w}{Vector of component weights.} \item{nodes}{List of nodes in
#' the subnetwork.} \item{K}{Number of mixture components.}
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references Leo Lahti et al.: Global modeling of transcriptional responses
#' in interaction networks. Bioinformatics (2010). See citation('netresponse')
#' for details.
#' @keywords utilities
#' @export
#' @examples
#' 
#' # Load toy data
#' data( toydata )          # Load toy data set
#' D     <- toydata$emat    # Response matrix (for example, gene expression)
#' model <- toydata$model   # Pre-calculated model
#' 
#' # Get model parameters for a given subnet
#' # (Gaussian mixture: mean, covariance diagonal, mixture proportions)
#' get.model.parameters(model, subnet.id = 1)
#' 
get.model.parameters <- function(model, subnet.id = NULL) {
    
    if (is(model) == "NetResponseModel") {
        
        if (is.numeric(subnet.id)) {
            warning("subnet.id given as numeric; converting to character: ", "Subnet-", 
                subnet.id, sep = "")
            subnet.id <- paste("Subnet", subnet.id, sep = "-")
        }
        
        pars <- model@models[[subnet.id]]
        pars[["nodes"]] <- model@subnets[[subnet.id]]
        
    } else {
        
        pars <- list()
        pars$mu <- model$model$posterior$centroids
        pars$sd <- model$model$posterior$sds
        pars$w <- model$model$posterior$weights
        
    }
    
    pars
    
}


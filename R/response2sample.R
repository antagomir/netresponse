#' response2sample
#' 
#' List the most strongly associated response of a given subnetwork for each
#' sample.
#' 
#' 
#' @param model A NetResponseModel object or list. 
#' @param subnet.id Subnet id. A natural number which specifies one of the
#' subnetworks within the 'model' object.
#' @param component.list List samples separately for each mixture component
#' (TRUE). Else list the most strongly associated component for each sample
#' (FALSE).
#' @param verbose Follow progress by intermediate messages.
#' @param data Data (features x samples; or a vector for univariate case) to predict response for given data points (currently implemented only for mixture.model output)
#'
#' Return:
#' @return A list. Each element corresponds to one subnetwork response, and
#' contains a list of samples that are associated with the response (samples
#' for which this response has the highest probability P(response | sample)).
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references Leo Lahti et al.: Global modeling of transcriptional responses
#' in interaction networks. Bioinformatics (2010).  See citation("netresponse")
#' for citation details.
#' @keywords utilities
#' @export
#' @examples
#' # Load example data
#' data( toydata )         # Load toy data set
#' D    <- toydata$emat    # Response matrix (for example, gene expression)
#' model <- toydata$model  # Pre-calculated model
#' 
#' # Find the samples for each response (for a given subnetwork)
#' response2sample(model, subnet.id = 1)
#' 
response2sample <- function (model, subnet.id = NULL, component.list = TRUE, verbose = FALSE, data = NULL) {

  if (class(model) == "NetResponseModel") {

    if (is.numeric(subnet.id)) {
      subnet.id <- paste("Subnet", subnet.id, sep = "-")
      warning("subnet.id given as numeric; converting to character: ", subnet.id, sep="")
    }
  
    response.probabilities <- model[[subnet.id]]$qofz # sample2response(model, subnet.id)

    rownames(response.probabilities) <- rownames(model@datamatrix)

  } else if (class(model) == "list") {

    if (!is.null(model$qofz)) {
      # Pick response probabilities from the model object
      response.probabilities <- model$qofz
    } else {
      # Otherwise, retrieve the response probabilities assuming the 
      # input data and parameters are presented on the same basis 
      # For mixture.model output
      if (is.vector(data)) { 
        data2 <- t(matrix(data)) 
        colnames(data2) <- names(data)
        rownames(data2) <- "Mode-1"
        data <- data2
        data2 <- NULL
      }

      # Find cluster for each sample
      response.probabilities <- P.r.s(data, model$params, log = TRUE)
      #rownames(response.probabilities) <- colnames(data)
      colnames(response.probabilities) <- paste("Mode", 1:ncol(response.probabilities), sep = "-")
    }
  }

  # For each sample, list the most strongly associated response (highest P(r|s))
  clusters <- apply(response.probabilities, 1, which.max)

  if (length(clusters) == 0) {
    warning(paste("Error with response2sample in subnet", subnet.id))
  } else if ( component.list ) {
    # list samples separately for each cluster
    clusters <- lapply(seq(max(clusters)), function( i ){ 
      if (is.null(names(clusters))) {
        which(clusters == i)
      } else {
        names(which(clusters == i))
      }
    })

    names(clusters) <- paste("Mode-", 1:length(clusters), sep = "")
    if ( length(clusters) < ncol(response.probabilities) ) {
      n <- ncol(response.probabilities) - length(clusters)
      clusters <- c(clusters, vector(n, mode = "list"))
    }
  }
    
  clusters
  
}


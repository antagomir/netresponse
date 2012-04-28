# Copyright (C) 2010-2012 Leo Lahti
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

#' response2sample
#' 
#' List the most strongly associated response of a given subnetwork for each
#' sample.
#' 
#' 
#' @usage response2sample(model, subnet.id, component.list = TRUE, verbose =
#' FALSE)
#' @param model A NetResponseModel object. Result from NetResponse
#' (detect.responses function).
#' @param subnet.id Subnet id. A natural number which specifies one of the
#' subnetworks within the 'model' object.
#' @param component.list List samples separately for each mixture component
#' (TRUE). Else list the most strongly associated component for each sample
#' (FALSE).
#' @param verbose Follow progress by intermediate messages.
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
#'  
#' 
#' library( netresponse )
#' 
#' # Load example data
#' data( toydata )         # Load toy data set
#' D    <- toydata$emat    # Response matrix (for example, gene expression)
#' model <- toydata$model  # Pre-calculated model
#' 
#' # Find the samples for each response (for a given subnetwork)
#' response2sample(model, subnet.id = 1)
#' 
#' 
response2sample <- function (model, subnet.id, component.list = TRUE, verbose = FALSE) {

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", subnet.id, sep="")
  }
  
  response.probabilities <- sample2response(model, subnet.id)

  # For each sample, list the most strongly associated response (highest P(r|s))
  clusters <- apply(response.probabilities, 1, which.max)

  if (length(clusters) == 0) {
    warning(paste("Error with response2sample in subnet", subnet.id))
  } else if ( component.list ) {
    # list samples separately for each cluster
    clusters <- lapply(seq(max(clusters)), function( i ){ names(which(clusters == i)) })
    # Names(clusters) <- FIXME add names here
    if ( length(clusters) < ncol(response.probabilities) ) {
      n <- ncol(response.probabilities) - length(clusters)
      clusters <- c(clusters, vector(n, mode = "list"))
    }
  }
    
  clusters
  
}


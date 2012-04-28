# Copyright (C) 2010-2012 Leo Lahti
# Contact: Leo Lahti <leo.lahti@iki.fi>
# This file is part of NetResponse R program
#  
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



#' sample2response
#' 
#' Probabilistic sample-response assignments for given subnet.
#' 
#' 
#' @usage sample2response(model, subnet.id, mode = "soft")
#' @param model Result from NetResponse (detect.responses function).
#' @param subnet.id Subnet identifier. A natural number which specifies one of
#' the subnetworks within the 'model' object.
#' @param mode soft: gives samples x responses probabilistic assignment matrix;
#' hard: gives the most likely response for each sample
#' @return A matrix of probabilities. Sample-response assignments for given
#' subnet, listing the probability of each response, given a sample.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references Leo Lahti et al.: Global modeling of transcriptional responses
#' in interaction networks. Bioinformatics (2010).  See citation("netresponse")
#' for citation details.
#' @keywords utilities
#' @export
#' @examples
#' 
#' 
#' #library(netresponse)
#' #data( toydata )        # Load toy data set
#' #D    <- toydata$emat   # Response matrix (for example, gene expression)
#' #netw <- toydata$netw   # Network
#' 
#' # Detect network responses
#' #model <- detect.responses(D, netw, verbose = FALSE)
#' 
#' # Assign samples to responses (soft, probabilistic assignments sum to 1)
#' #response.probabilities <- sample2response(model, subnet.id = "Subnet-1")
#' 
#' 
sample2response <- function (model, subnet.id, mode = "soft") {

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", subnet.id, sep="")
  }
  
  # P(response | sample)
  #assignment.matrix <- model@models[[subnet.id]]$posterior$qOFz
  assignment.matrix <- getqofz(model, subnet.id, log = FALSE)
  
  if (mode == "hard") {
    sample.names <- rownames(assignment.matrix)
    assignment.matrix <- colnames(assignment.matrix)[apply(assignment.matrix, 1, which.max)]
    names(assignment.matrix) <- sample.names
  }

  assignment.matrix

}



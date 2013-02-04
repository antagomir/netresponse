# Copyright (C) 2010-2013 Leo Lahti
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


#' model.stats
#' 
#' Subnetwork statistics: size and number of distinct responses for each
#' subnet.
#' 
#' 
#' @param models NetResponse object or list of models
#' @return
#' 
#' A 'subnetworks x properties' data frame containing the following elements.
#' \item{subnet.size: }{ Vector of subnetwork sizes. } \item{subnet.responses:
#' }{ Vector giving the number of responses in each subnetwork. }
#' @author Leo Lahti <leo.lahti@@iki.fi>
#' @references Leo Lahti et al.: Global modeling of transcriptional responses
#' in interaction networks. Bioinformatics (2010).  See citation("netresponse")
#' for reference details.
#' @keywords utilities
#' @export
#' @examples
#' 
#' 
#' library(netresponse)
#' 
#' # Load a pre-calculated netresponse model obtained with 
#' # model <- detect.responses(toydata$emat, toydata$netw, verbose = FALSE)
#' data( toydata )        
#' # Calculate summary statistics for the model
#' stat <- model.stats(toydata$model)
#' 
#' 
model.stats <- function ( models ) {

  # Check statistics for subnetworks
  # subnetwork size
  # number of responses

  if (class(models) == "NetResponseModel") { 
    models <- models@models
  }

  if (is.null(names(models))) {names(models) <- 1:length(models)}

  Ncomps <- c()
  for (subnet.id in names(models)) {
    Ncomps[[subnet.id]] <- ncol(models[[subnet.id]]$qofz)
  }

  tab <- cbind(sapply(models, function (x) {ncol(x$mu)}), Ncomps)
  colnames(tab) <- c("nodes", "responses")
  rownames(tab) <- names(models)

  as.data.frame(tab)

}


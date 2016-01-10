# Copyright (C) 2010-2016 Leo Lahti
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


# "To invent, you need a good imagination and a pile of junk." 
#     	      	       	      -- Thomas Edison

#' @title Continuous responses
#' @description Quantify association between modes and continuous variable
#' @param annotation.vector annotation vector with discrete factor levels, and named by the samples
#' @param model NetResponse model object
#' @param method method for enrichment calculation
#' @param min.size minimum sample size for a response 
#' @param data data matrix (samples x features)
#' @return List with each element corresponding to one variable and listing the responses according to association strength
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
continuous.responses <- function (annotation.vector, model, method = "t-test", min.size = 2, data = NULL) {
		    
  # method = "t-test"; min.size = 2; data = t(dat[, gpt])		     

  if (is.null(data) && class(model) == "NetResponseModel") {
    data <- model@datamatrix
    all.samples <- rownames(data)
  } 

  # samples x features
  if(is.vector(data)) {
    data2 <- matrix(data)
    rownames(data2) <- names(data)
    data <- data2
    all.samples <- rownames(data)
  }

  if (class(model) == "NetResponseModel") {
    models <- get.subnets(model, min.size = min.size)
  } else {
    models <- model
  }

  associations <- enrichment.list(models, annotation.vector)

  associations

}

#associations <- continuous.responses.single(model, annotation.vector, method = "t.test")
continuous.responses.single <- function (model, annotation.vector, method = "t.test") {

  annotated.samples <- names(which(!is.na(annotation.vector)))
  annotation.data <- annotation.vector[annotated.samples]
  names(annotation.data) <- annotated.samples


  r2s <- response2sample(model)

  if (length(r2s) <= 1) {return(NA)} # No multimodality -> no enrichments

  pvals <- c()
  fold.change <- c()

  for (mo in 1:length(r2s)) {

    # annotated samples in the mode
    s <- intersect(r2s[[mo]], annotated.samples)

    # annotated samples in other modes
    sc <- intersect(unlist(r2s[-mo]), annotated.samples)

    if (length(na.omit(s)) > 1 && length(na.omit(sc)) > 1) {     
      if (method == "t.test") {
        pval <- t.test(annotation.data[s], annotation.data[sc])$p.value 
      }

      pvals[[mo]] <- pval

      fold.change[[mo]] <- mean(annotation.data[s]) - mean(annotation.data[sc])
    } else {

      warning(paste("Not enough annotated observations to calculate p-values", mo))
      pvals[[mo]] <- NA
      fold.change[[mo]] <- NA
    }
  }

  associations <- data.frame(list(mode = 1:length(r2s), pvalue = pvals, fold.change = fold.change))

  associations
      
}





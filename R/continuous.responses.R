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


# "To invent, you need a good imagination and a pile of junk." 
#     	      	       	      -- Thomas Edison


#' Description: Quantify association between modes and continuous variable
#' 
#' Arguments:
#'   @param annotation.vector annotation vector with discrete factor levels, and named by the samples
#'   @param model NetResponse model object
#'   @param method method for enrichment calculation
#'   @param min.size minimum sample size for a response 
#'   @param data data matrix (samples x features)
#'
#' Returns:
#'   @return List with each element corresponding to one variable and listing the responses according to association strength
#'            
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
continuous.responses <- function (annotation.vector, model, method = "t-test", min.size = 2, data = NULL) {
		    
  # method = "t-test"; min.size = 2; data = t(dat[, gpt])		     

  if (is.null(data) && class(model) == "NetResponseModel") {
    data <- model@datamatrix
  } else if (is.null(data)) {
    stop("Provide data")
  }

  # samples x features
  if(is.vector(data)) {
    data2 <- matrix(data)
    rownames(data2) <- names(data)
    data <- data2
  }

  all.samples <- rownames(data)
  annotated.samples <- names(which(!is.na(annotation.vector)))
  annotation.data <- annotation.vector[annotated.samples]
  names(annotation.data) <- annotated.samples

  associations <- NULL  	 
  if (class(model) == "NetResponseModel") {

    subnets <- get.subnets(model, min.size = min.size)

    for (sn in names(subnets)) {

      # samples in each mode (hard assignment)
      r2s <- response2sample(model, subnet.id = sn)

      pvals <- c()
      mean.difference <- c()
      for (mo in 1:length(r2s)) {

        # annotated samples in the mode
        s <- intersect(r2s[[mo]], annotated.samples)

	# annotated samples in other modes
        sc <- intersect(unlist(r2s[-mo]), annotated.samples)

        if (length(na.omit(s)) > 1 && length(na.omit(sc)) > 1) {      
          pvals[[mo]] <- t.test(annotation.data[s], annotation.data[sc])$p.value
          mean.difference[[mo]] <- mean(annotation.data[s]) - mean(annotation.data[sc])

        } else {
          warning(paste("Not enough annotated observations for response", mo))
          pvals[[mo]] <- NA
	  mean.difference[[mo]] <- NA
        }
      }

      associations <- rbind(associations, cbind(subnet = rep(sn, length(r2s)), mode = 1:length(r2s), pvalue = pvals, mean.difference = mean.difference))

    }

    associations <- data.frame(list(subnet = associations[, "subnet"], 
    		    	            mode = associations[, "mode"], 
    		    	            pvalue = as.numeric(associations[, "pvalue"]), 
    		    	            mean.difference = as.numeric(associations[, "mean.difference"])
    				    ))

  } else if (class(model) == "list") {

    # for mixture.model output, for instance; assuming there is only a single 'subnet'
 
    if (all(c("mu", "sd", "w") %in% names(model))) {
      # Assuming it is a single group
      associations <- continuous.responses.single(model, data, annotation.data, annotated.samples, method = "t.test")

    } 

  }

  associations

}


continuous.responses.single <- function (model, data, annotation.data, annotated.samples, method = "t.test") {

  r2s <- response2sample(model)

  pvals <- c()
  mean.difference <- c()

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

      mean.difference[[mo]] <- mean(annotation.data[s]) - mean(annotation.data[sc])
    } else {

      warning(paste("Not enough annotated observations to calculate p-values", mo))
      pvals[[mo]] <- NA
      mean.difference[[mo]] <- NA
    }
  }

  associations <- data.frame(list(mode = 1:length(r2s), pvalue = pvals, mean.difference = mean.difference))

  associations
      
}
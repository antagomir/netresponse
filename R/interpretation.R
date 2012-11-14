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


# "To invent, you need a good imagination and a pile of junk." 
#     	      	       	      -- Thomas Edison


#' Description: List responses for each level of the given factor
#' 
#' Arguments:
#'   @param annotation.vector annotation vector with discrete factor levels, and named by the samples
#'   @param model NetResponse model object
#'   @param method method for enrichment calculation
#'   @param min.size minimum sample size for a response 
#'   @param data data (samples x features; or a vector in univariate case)
#'
#' Returns:
#'   @return List with each element corresponding to one factor level and listing the responses according to association strength
#'            
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
factor.responses <- function (annotation.vector, model, method = "hypergeometric", min.size = 2, data = NULL) {

  responses <- list()

  annotation.vector <- factor(annotation.vector)

  levels <- as.character(na.omit(unique(droplevels(annotation.vector))))

  for (lev in levels) {

    level.samples <- names(annotation.vector)[which(annotation.vector == lev)]

      ors <- order.responses(model, level.samples, method = method, min.size = min.size, data = data) 
 
      if (is.null(ors)) { 
        ors <- NA 
        warning(paste("No significant responses for level", lev))
      }
  
    responses[[as.character(lev)]] <- ors
  }

  # Pick top responses for each factor level
  responses <- responses[!is.na(responses)]

  responses.per.level <- NULL
  if (length(responses) > 0) {

    responses.per.level <- lapply(responses, function (dr) { dr$ordered.responses })
    responses.per.level <- responses.per.level[sapply(responses.per.level, nrow) > 0]

  }

  responses.per.level

}






#' Description: List responses for all factors and levels in the given
#' annotation matrix
#' 
#' Arguments:  
#' @param annotation.df annotation data.frame with discrete factor levels, rows
#'   named by the samples
#' @param model NetResponse model object
#' @param method method for enrichment calculation
#' @param min.size minimum sample size for a response
#' @param qth q-value threshold
#' @param verbose verbose 
#' @param data data (samples x features; or a vector in univariate case) 
#' @param rounding rounding digits 
#'
#' Returns:
#' @return Table listing all associations between the factor levels and
#'   responses
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities

list.responses.factor <- function (annotation.df, model, method = "hypergeometric", min.size = 2, qth = Inf, verbose = TRUE, data = NULL, rounding = NULL) {

  # annotation.df <- atlas.metadata[sample.set, factor.vars]; model <- res$model; method = "hypergeometric"; min.size = 1; qth = Inf; verbose = TRUE
  # annotation.df <- annot[, factor.vars]; model; min.size = 1; qth = 1; method = "hypergeometric"; verbose = TRUE
  # annotation.df <- annot[sample.set, factor.vars]; model; min.size = 1; qth = qth; data = t(X); method = "hypergeometric"; verbose = TRUE

  # samples x features
  if(is.vector(data)) {
    data2 <- matrix(data)
    rownames(data2) <- names(data)
    data <- data2
  }

  if (!is.data.frame(annotation.df)) {
    stop("Provide data.frame for the annotation.df argument!")
  }

  # Collect the tables from all factors and levels here
  collected.table <- NULL

  for (fnam in colnames(annotation.df)) { 

    if (verbose) { message(fnam) }

    annotation.vector <- annotation.df[[fnam]]
    names(annotation.vector) <- rownames(annotation.df)

    # Order responses for each level of this factor
    responses.per.level <- factor.responses(annotation.vector, model, method = method, min.size = min.size, data = data) 
    responses.per.level <- responses.per.level[na.omit(names(responses.per.level))]

    # Add factor/level information
    mat <- NULL
    if (length(responses.per.level) > 0 && !is.null(na.omit(names(responses.per.level)))) {
      for (level in na.omit(names(responses.per.level))) {

        responses.per.level[[level]] <- cbind(responses.per.level[[level]], 
      				   rep(fnam,  nrow(responses.per.level[[level]])), 
				   rep(level, nrow(responses.per.level[[level]])))

        tmp <- responses.per.level[[level]]
        colnames(responses.per.level[[level]]) <- c(colnames(tmp)[1:(ncol(tmp)-2)], "Factor", "Level")

      }
    

      # Combine the level-wise matrices
      mat <- responses.per.level[[1]]; 

      if (length(responses.per.level) > 1) {
        for (i in 2:length(responses.per.level)) { mat <- rbind(mat, responses.per.level[[i]])} 
      }

    }

    # FIXME: do.call(rbind, l) tai do.call(cbind, l)
    # is a neater way to implement this
    collected.table <- rbind(collected.table, mat)

  }

  # Order by pvalue
  collected.table <- collected.table[order(as.numeric(as.character(collected.table[, "pvalue"]))), ]
  collected.table <- data.frame(collected.table)
  collected.table$pvalue <- as.numeric(as.character(collected.table$pvalue))

  if (!is.null(collected.table)) {
  
    collected.table$qvalue <- qvalue::qvalue(collected.table$pvalue, gui = FALSE)$qvalue
    colnames(collected.table) <- c(colnames(collected.table)[1:(ncol(collected.table)-1)], "qvalue")
  
    # Filtering based on qvalues
    collected.table <- collected.table[collected.table$qvalue < qth, ]
  
  }


  collected.table$mode <- as.character(collected.table$mode)
  #collected.table$Factor <- collected.table$Factor
  #collected.table$Level <- collected.table$Level
  collected.table$mysamples.in.response <- as.numeric(as.character(collected.table$mysamples.in.response))
  collected.table$fraction.in.response <- as.numeric(as.character(collected.table$fraction.in.response))
  collected.table$fraction.in.data <- as.numeric(as.character(collected.table$fraction.in.data))
  collected.table$pvalue <- as.numeric(as.character(collected.table$pvalue))
  collected.table$qvalue <- as.numeric(as.character(collected.table$qvalue))
 

  if (!is.null(rounding)) {

    collected.table$fraction.in.response <- round(collected.table$fraction.in.response, rounding)
    collected.table$fraction.in.data <- round(collected.table$fraction.in.data, rounding)
    collected.table$pvalue <- round(collected.table$pvalue, rounding)
    collected.table$qvalue <- round(collected.table$qvalue, rounding)
 
  }
 
  collected.table

}

#################################################################################################


#' Description: Investigate association of a continuous variable and the modes
#' 
#' Arguments: 
#' @param annotation.df annotation data.frame with discrete factor levels, rows
#' named by the samples
#' @param model NetResponse model object
#' @param method method for quantifying the association
#' @param min.size minimum sample size for a response
#' @param qth q-value threshold
#' @param verbose verbose 
#' @param data data (samples x features)
#' @param rounding rounding digits
#'
#' Returns:
#' @return Table listing all associations between the factor levels and
#' responses
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities

list.responses.continuous <- function (annotation.df, model, method = "t-test", min.size = 1, qth = Inf, verbose = TRUE, data = NULL, rounding = NULL) {

  # annotation.df <- atlas.metadata[sample.set, continuous.vars]; method <- "t-test"; model <- res$model; min.size = 1; qth = 0.2; verbose = TRUE
  # annotation.df <- annot[, continuous.vars];  method <- "t-test"; min.size = 1; qth = 0.2; verbose = TRUE
  # annotation.df <- annot[, continuous.vars]; method = "t-test"; min.size = 1; qth = qth; verbose = TRUE
  # source("~/Rpackages/netresponse/netresponse/R/interpretation.R")
  # annotation.df <- annot[, continuous.vars]; method = "t-test"; model; min.size = 1; qth = 1; data = z; verbose = TRUE

  # Collect the tables from all factors and levels here
  collected.table <- NULL

  for (fnam in colnames(annotation.df)) { 

    if (verbose) { message(fnam) }

    annotation.vector <- annotation.df[[fnam]]
    names(annotation.vector) <- rownames(annotation.df)

    # quantify association to each response for the continuous variable
    responses.per.cont <- continuous.responses(annotation.vector, model, method = method, min.size = min.size, data = data)
    responses.per.cont$annotation <- rep(fnam, nrow(responses.per.cont))

    collected.table <- rbind(collected.table, responses.per.cont)

  }

  if (nrow(collected.table) > 0) {

    collected.table$qvalue <- rep(NA, nrow(collected.table))
    nainds <- is.na(collected.table$pvalue)
    if (sum(!nainds) > 50) {
      qv <- qvalue(collected.table$pvalue[!nainds], pi0.method = "bootstrap")
      if (("qvalues" %in% names(qv)) && sum(!nainds) > 0) {
        collected.table$qvalue[!nainds] <- qv$qvalues
      } else {
        warning("Too few values for qvalue estimation, skipped")
        collected.table$qvalue[!nainds] <- rep(NA, sum(!nainds))
      }
    } 

    if (sum(nainds) > 0) {
      warning("Removing entries where p/q values could not be calculated due to small sample size and/or missing values")
      collected.table <- collected.table[!nainds,]
    }

    # Order by pvalues
    collected.table <- collected.table[order(collected.table$pvalue),]

    # Filtering based on qvalues, if qvalues are available
    if (any(!is.na(collected.table$qvalue)) && !is.null(qth)) {
      collected.table <- collected.table[collected.table$qvalue < qth, ]
    } 
  } else {
    collected.table <- NULL
  }

  if (length(collected.table) == 0) { collected.table <- NULL} 


  if (!is.null(rounding)) {
    collected.table$qvalue <- round(collected.table$qvalue, rounding)
    collected.table$pvalue <- round(collected.table$pvalue, rounding)
  }  


  collected.table

}



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

  # annotation.vector, model, method = method, min.size = min.size

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
      for (mo in 1:length(r2s)) {

        # annotated samples in the mode
        s <- intersect(r2s[[mo]], annotated.samples)

	# annotated samples in other modes
        sc <- intersect(unlist(r2s[-mo]), annotated.samples)

        if (length(na.omit(s)) > 1 && length(na.omit(sc)) > 1) {      
          pvals[[mo]] <- t.test(annotation.data[s], annotation.data[sc])$p.value
        } else {
          warning(paste("Not enough annotated observations for response", mo))
          pvals[[mo]] <- NA
        }
    }

    associations <- rbind(associations, cbind(subnet = rep(sn, length(r2s)), mode = paste("Mode-", 1:length(r2s), sep = ""), pvalue = pvals))

    }

  } else if (class(model) == "list") {

    # for mixture.model output, for instance; assuming there is only a single 'subnet'
         
    # samples in each mode (hard assignment)
    #r2s <- model$model$posterior$qOFz 
    r2s <- response2sample(model, data = t(data))

    pvals <- c()
    for (mo in 1:length(r2s)) {

      # annotated samples in the mode
      s <- intersect(r2s[[mo]], annotated.samples)

      # annotated samples in other modes
      sc <- intersect(unlist(r2s[-mo]), annotated.samples)

      if (length(na.omit(s)) > 1 && length(na.omit(sc)) > 1) {      
        pvals[[mo]] <- t.test(annotation.data[s], annotation.data[sc])$p.value
      } else {
        warning(paste("Not enough annotated observations for response", mo))
        pvals[[mo]] <- NA
      }
    }

    associations <- data.frame(list(mode = paste("Mode-", 1:length(r2s), sep = ""), pvalue = pvals))

  }

  associations

}

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
#'
#' Returns:
#'   @return List with each element corresponding to one factor level and listing the responses according to association strength
#'            
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
factor.responses <- function (annotation.vector, model, method = "hypergeometric", min.size = 2) {

  # annotation.vector, model, method = method, min.size = min.size

  responses <- list()

  annotation.vector <- factor(annotation.vector)

  levels <- as.character(na.omit(unique(droplevels(annotation.vector))))

  for (lev in levels) {
    level.samples <- names(annotation.vector)[which(annotation.vector == lev)]
    ors <- order.responses(model, level.samples, method = method, min.size = min.size) 
    if (is.null(ors)) { 
      ors <- NA 
      warning(paste("No significant responses for level", lev))
    }
    responses[[as.character(lev)]] <- ors
  }

  # Pick top responses for each factor level

  responses.per.level <- lapply(responses, function (dr) { dr$ordered.responses })
  responses.per.level <- responses.per.level[sapply(responses.per.level, nrow) > 0]

  responses.per.level

}




#' Description: Quantify association between modes and continuous variable
#' 
#' Arguments:
#'   @param annotation.vector annotation vector with discrete factor levels, and named by the samples
#'   @param model NetResponse model object
#'   @param method method for enrichment calculation
#'   @param min.size minimum sample size for a response 
#'
#' Returns:
#'   @return List with each element corresponding to one factor level and listing the responses according to association strength
#'            
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
continuous.responses <- function (annotation.vector, model, method = "t-test", min.size = 2) {

  # annotation.vector, model, method = method, min.size = min.size

  subnets <- get.subnets(model, min.size = min.size)

  associations <- NULL  	 

  for (sn in names(subnets)) {

    # samples in each mode (hard assignment)
    r2s <- response2sample(model, subnet.id = sn)

    all.samples <- rownames(model@datamatrix)
    ann.inds <- !is.na(annotation.vector)
    annotation.samples <- all.samples[ann.inds]
    annotation.data <- annotation.vector[ann.inds]
    names(annotation.data) <- annotation.samples

    pvals <- c()
    for (mo in 1:length(r2s)) {
      # Only consider annotated samples
      s <- intersect(r2s[[mo]], annotation.samples)
      sc <- setdiff(annotation.samples, s)
      if (length(s) > 1 && length(sc) > 1) {      
        pvals[[mo]] <- t.test(annotation.data[s], annotation.data[sc])$p.value
      } else {
         warning(paste("Not enough annotated observations for response", mo))
        pvals[[mo]] <- NA
      }
    }

    associations <- rbind(associations, cbind(subnet = rep(sn, length(r2s)), mode = paste("Mode-", 1:length(r2s), sep = ""), pvalue = pvals))

  }

  associations <- data.frame(list(subnet = associations[, "subnet"], mode = associations[, "mode"], pvalue = as.numeric(associations[, "pvalue"])))

  nainds <- is.na(associations$pvalue)

  associations$qvalue <- rep(NA, nrow(associations))
  if (sum(!nainds) > 20) {
    associations$qvalue[!nainds] <- qvalue(associations$pvalue[!nainds])$qvalue
  } else {
    warning("Not enough pvalues for qvalue estimation, skipping.")
  }

  associations <- associations[order(associations$qvalue), ]

  associations

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
#' @param verbose verbose Returns:
#' @return Table listing all associations between the factor levels and
#'   responses
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities

list.responses.factor <- function (annotation.df, model, method = "hypergeometric", min.size = 2, qth = Inf, verbose = TRUE) {

  # annotation.df <- atlas.metadata[sample.set, factor.vars]; model <- res$model; method = "hypergeometric"; min.size = 1; qth = Inf; verbose = TRUE

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
    responses.per.level <- factor.responses(annotation.vector, model, method = method, min.size = min.size) 
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
        responses.per.level[[level]][, "qvalue"] <- NULL
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

  if (!is.null(collected.table)) {

    collected.table <- cbind(collected.table, qvalue::qvalue(collected.table[, "pvalue"], gui = FALSE)$qvalue)
    colnames(collected.table) <- c(colnames(collected.table)[1:(ncol(collected.table)-1)], "qvalue")

    # Order by qvalue
    collected.table <- collected.table[order(collected.table$qvalue), ]

    # Filtering based on qvalues
    collected.table[collected.table$qvalue < qth, ]

  }
  
}



#' Description: Investigate association of a continuous variable and the modes
#' 
#' Arguments: 
#' @param annotation.df annotation data.frame with discrete factor levels, rows
#' named by the samples
#' @param model NetResponse model object
#' @param method method for quantifying the association
#' @param min.size minimum sample size for a response
#' @param qth q-value threshold
#' @param verbose verbose Returns:
#' @return Table listing all associations between the factor levels and
#' responses
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities

list.responses.continuous <- function (annotation.df, model, method = "t-test", min.size = 1, qth = Inf, verbose = TRUE) {

  # annotation.df <- annot[, continuous.vars]; method = "t-test"; min.size = 1; qth = qth; verbose = TRUE

  # Collect the tables from all factors and levels here
  collected.table <- NULL

  for (fnam in colnames(annotation.df)) { 

    if (verbose) { message(fnam) }

    annotation.vector <- annotation.df[[fnam]]
    names(annotation.vector) <- rownames(annotation.df)

    # quantify association to each response for the continuous variable
    responses.per.cont <- continuous.responses(annotation.vector, model, method = method, min.size = min.size)
    responses.per.cont$annotation <- rep(fnam, nrow(responses.per.cont))

    collected.table <- rbind(collected.table, responses.per.cont)

  }

  collected.table$qvalue <- rep(NA, nrow(collected.table))
  nainds <- is.na(collected.table$pvalue)
  collected.table$qvalue[!nainds] <- qvalue(collected.table$pvalue[!nainds])$qvalue
  
  if (sum(nainds) > 0) {
    warning("Removing entries where p/q values could not be calculated due to small sample size and/or missing values")
    collected.table <- collected.table[!nainds,]
  }

  # Order by qvalues
  collected.table <- collected.table[order(collected.table$qvalue),]

  # Filtering based on qvalues
  collected.table[collected.table$qvalue < qth, ]


}



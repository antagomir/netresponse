<<<<<<< HEAD
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

#' Description: List significantly associated responses for all factors and levels in the given
#' annotation matrix
#' 
#' Arguments:  
#' @param annotation.df annotation data.frame with discrete factor levels, rows
#'   named by the samples
=======
#' @title List significant responses
#' @description List significantly associated responses for all factors and levels in the given annotation matrix
#' @param annotation.df annotation data.frame with discrete factor levels, rows named by the samples
>>>>>>> master
#' @param models List of models. Each model should have a sample-cluster assignment matrix qofz, or a vector of cluster indices named by the samples.
#' @param method method for enrichment calculation
#' @param min.size minimum sample size for a response
#' @param qth q-value threshold
#' @param verbose verbose 
#' @param data data (samples x features; or a vector in univariate case) 
#' @param rounding rounding digits 
<<<<<<< HEAD
#'
#' Returns:
#' @return Table listing all associations between the factor levels and
#'   responses
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities

=======
#' @return Table listing all associations between the factor levels and responses
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @importFrom qvalue qvalue
#' @export
#' @keywords utilities
>>>>>>> master
list.responses.factor <- function (annotation.df, models, method = "hypergeometric", min.size = 2, qth = Inf, verbose = TRUE, data = NULL, rounding = NULL) {

  pth <- NULL

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
    responses.per.level <- factor.responses(annotation.vector, models, method = method, min.size = min.size, data = data) 

    responses.per.level <- responses.per.level[na.omit(names(responses.per.level))]

    # Add factor/level information
    mat <- NULL
    if (length(responses.per.level) > 0 && !is.null(na.omit(names(responses.per.level)))) {
      for (level in na.omit(names(responses.per.level))) {

        responses.per.level[[level]] <- cbind(responses.per.level[[level]], 
							fnam, 
				   			level)
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
  
    if (nrow(collected.table)>100) {
<<<<<<< HEAD
      collected.table$qvalue <- qvalue::qvalue(collected.table$pvalue, gui = FALSE, fdr.level = 0.25)$qvalue
=======
      collected.table$qvalue <- qvalue(collected.table$pvalue, gui = FALSE, fdr.level = 0.25)$qvalue
>>>>>>> master
    } else {
      collected.table$qvalue <- rep(NA, nrow(collected.table))
    }  

    # Filtering based on qvalues
    if (!all(is.na(collected.table$qvalue))) {
      collected.table <- collected.table[collected.table$qvalue < qth, ]
    } else {
      collected.table <- collected.table[collected.table$pvalue < pth, ]
    }
  }

  collected.table$mode <- as.character(collected.table$mode)
  collected.table$annotated.in.subset <- as.numeric(as.character(collected.table$annotated.in.subset))
  collected.table$fraction.in.subset <- as.numeric(as.character(collected.table$fraction.in.subset))
  collected.table$fraction.in.data <- as.numeric(as.character(collected.table$fraction.in.data))
  collected.table$pvalue <- as.numeric(as.character(collected.table$pvalue))
  collected.table$qvalue <- as.numeric(as.character(collected.table$qvalue))
 
  if (!is.null(rounding)) {

    collected.table$fraction.in.subset <- round(collected.table$fraction.in.subset, rounding)
    collected.table$fraction.in.data <- round(collected.table$fraction.in.data, rounding)
    collected.table$pvalue <- round(collected.table$pvalue, rounding)
    collected.table$qvalue <- round(collected.table$qvalue, rounding)
 
  }

  collected.table

}


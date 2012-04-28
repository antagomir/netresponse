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

#' Description: List responses for each level of the given factor
#' 
#' Arguments:
#'   @param annotation.vector annotation vector with discrete factor levels, and named by the samples
#'   @param model NetResponse model object
#'   @param method method for enrichment calculation
#'   @param min.size minimum sample size for a response 
#'   @param qth q-value threshold
#'
#' Returns:
#'   @return List with each element corresponding to one factor level and listing the responses according to association strength
#'            
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
factor.responses <- function (annotation.vector, model, method = "hypergeometric", min.size = 2, qth = Inf) {

  responses <- list()
  levels <- unique(annotation.vector)
  for (lev in levels) {
    sample <- names(annotation.vector)[annotation.vector == lev]
    responses[[lev]] <- order.responses(model, sample, method = method, min.size = min.size) 
  }

  # Pick top responses for each factor level
  responses.per.level <- lapply(responses, function (dr) { subset(dr$ordered.responses, qvalue < qth) })
  responses.per.level <- responses.per.level[sapply(responses.per.level, nrow) > 0]

  responses.per.level

}


#' Description: List responses for all factors and levels in the given
#' annotation matrix
#' 
#' Arguments:
#' 
#' 
#' @usage list.responses(annotation.df, model, method = "hypergeometric",
#' min.size = 2, qth = Inf, verbose = TRUE)
#' @param annotation.df annotation data.frame with discrete factor levels, rows
#' named by the samples
#' @param model NetResponse model object
#' @param method method for enrichment calculation
#' @param min.size minimum sample size for a response
#' @param qth q-value threshold
#' @param verbose verbose Returns:
#' @return Table listing all associations between the factor levels and
#' responses
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
list.responses <- function (annotation.df, model, method = "hypergeometric", min.size = 2, qth = Inf, verbose = TRUE) {

  # Collect the tables from all factors and levels here
  collected.table <- NULL

  for (fnam in colnames(annotation.df)) { 

    if (verbose) { message(fnam) }

    annotation.vector <- annotation.df[[fnam]]
    names(annotation.vector) <- rownames(annotation.df)

    # Order responses for each level of this factor
    responses.per.level <- factor.responses(annotation.vector, model, method = method, min.size = min.size, qth = Inf) 
    responses.per.level <- responses.per.level[na.omit(names(responses.per.level))]

    # Add factor/level information
    for (level in na.omit(names(responses.per.level))) {
      
      responses.per.level[[level]] <- cbind(responses.per.level[[level]], 
      				   rep(fnam, nrow(responses.per.level[[level]])), 
				   rep(level, nrow(responses.per.level[[level]])))

      tmp <- responses.per.level[[level]]
      colnames(responses.per.level[[level]]) <- c(colnames(tmp)[1:(ncol(tmp)-2)], "Factor", "Level")
      responses.per.level[[level]][, "qvalue"] <- NULL

    }

    # Combine the level-wise matrices
    mat <- responses.per.level[[1]]; 
    for (i in 2:length(responses.per.level)) { mat <- rbind(mat, responses.per.level[[i]])} 
    # FIXME: do.call(rbind, l) tai do.call(cbind, l)
    # is a neater way to implement this

    collected.table <- rbind(collected.table, mat)

  }

  collected.table <- cbind(collected.table, qvalue(collected.table[, "pvalue"])$qvalue)
  colnames(collected.table) <- c(colnames(collected.table)[1:(ncol(collected.table)-1)], "qvalue")
  
  collected.table

}




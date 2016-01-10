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

#' @title Investigate association of a continuous variable and the modes
#' @description Investigate association of a continuous variable and the modes.
#' @param annotation.df annotation data.frame with discrete factor levels, rows
#' named by the samples
#' @param groupings Sample mode information. Each element corresponds to one of the modes and lists the samples assignment matrix qofz. Alternatively, a vector of mode indices named by the samples can be given.
#' @param method method for quantifying the association
#' @param pth p-value threshold (for adjusted p-values)
#' @param verbose verbose 
#' @param rounding rounding digits
#' @param adjust.p Adjust p-values (this will add p.adj column and remove pvalue column in the output table)
#' @return Table listing all associations between the factor levels and responses
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @importFrom qvalue qvalue
#' @export
#' @keywords utilities
list.responses.continuous.single <- function (annotation.df, groupings, method = "t-test", pth = Inf, verbose = TRUE, rounding = NULL, adjust.p = TRUE) {

  # Collect the tables from all factors and levels here
  collected.table <- NULL

  # Quantify association to each response for the continuous variable
  if (is.null(names(groupings))) {names(groupings) <- 1:length(groupings)}

  if (verbose) {message("Check mode lists")}

  # Convert grouping info to a list
  groupings.list <- listify.groupings(groupings)

  if (verbose) {message("Go through annotations")}
  for (fnam in colnames(annotation.df)) { 

    if (verbose) { message(fnam) }

    annotation.vector <- annotation.df[[fnam]]
    names(annotation.vector) <- rownames(annotation.df)

    responses.per.cont <- enrichment.list(groupings.list, annotation.vector)

    responses.per.cont$annotation <- rep(fnam, nrow(responses.per.cont))

    collected.table <- rbind(collected.table, responses.per.cont)

  }

  if (nrow(collected.table) > 0) {

    nainds <- is.na(collected.table$pvalue)

    if (sum(nainds) > 0) {
      warning("Removing entries where p-values could not be calculated due to small sample size and/or missing values")
      collected.table <- collected.table[!nainds,]
    }

    if (adjust.p) {

      collected.table$p.adj <- rep(NA, nrow(collected.table))
      if (nrow(collected.table) > 100) {
        if (verbose) {message("Adjusting p with q")}
        qv <- qvalue(collected.table$pvalue, pi0.method = "bootstrap", fdr.level = 0.25)
        if (("qvalues" %in% names(qv))) {
          collected.table$p.adj <- qv$qvalues
        }
      } else {
        if (verbose) {message("Adjusting p with BH")}
        collected.table$p.adj <- p.adjust(collected.table$pvalue, method = "BH")
      }
    }

    # Order by pvalues
    collected.table <- collected.table[order(collected.table$pvalue),]

    # Filtering based on p.adjs, if p.adjs are available
    if (adjust.p && (any(!is.na(collected.table$p.adj)) && !is.null(pth))) {
      collected.table <- collected.table[collected.table$p.adj < pth, ]
    } 
  } else {
    collected.table <- NULL
  }
  if (length(collected.table) == 0) { collected.table <- NULL} 

  if (!is.null(rounding)) {
    collected.table$p.adj <- round(collected.table$p.adj, rounding)
    collected.table$pvalue <- round(collected.table$pvalue, rounding)
    collected.table$fold.change <- round(collected.table$fold.change, rounding)
  }  

  collected.table

}


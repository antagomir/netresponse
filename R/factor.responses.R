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

#' @title Factor responses
#' @description List responses for each level of the given factor
#' @param annotation.vector annotation vector with discrete factor levels, and named by the samples
#' @param groupings List of groupings. Each model should have a sample-cluster assignment matrix qofz, or a vector of cluster indices named by the samples.
#' @param method method for enrichment calculation
#' @param min.size minimum sample size for a response 
#' @param data data (samples x features; or a vector in univariate case)
#' @return List with each element corresponding to one factor level and listing the responses according to association strength
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
factor.responses <- function (annotation.vector, groupings, method = "hypergeometric", min.size = 2, data = NULL) {

  responses <- list()

  annotation.vector <- factor(annotation.vector)

  levels <- as.character(na.omit(unique(droplevels(annotation.vector))))

  for (lev in levels) {

    level.samples <- names(annotation.vector)[which(annotation.vector == lev)]

    ors <- enrichment.list.factor(groupings, level.samples, method = method)

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



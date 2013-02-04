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

#' Enrichment for a specified sample group in the given response.
#' 
#' Calculate enrichment values for a specified sample group in the given
#' response.
#'
#' @param qofz sample-mode assignment matrix
#' @param annotation.sample User-defined sample group. For instance, samples belonging to a particular annotation class.
#' @param response Response id (integer) within the subnet.
#' @param method Enrichment method.
#'
#' @return List with enrichment statistics, depending on enrichment method.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso order.responses
#' @references See citation("netresponse").
#' @keywords utilities
#' @export
#' @examples #enr <- response.enrichment(subnet.id, models, sample, response, method)
#' 
response.enrichment <- function (qofz, annotation.sample, response, method = "hypergeometric") {

  # pick sample data for the response and
  # ensure this is a matrix also when a single sample is given
  if (any(!annotation.sample %in% rownames(qofz))) {
    warning("Not all annotation samples are in the original data matrix; only the shared ones used for enrichment analysis.")
    annotation.sample <- intersect(annotation.sample, rownames(qofz))
  }

  response.samples <- response2sample(list(qofz = qofz), component.list = TRUE)
  if (length(response.samples) == 0) { warning("No samples in response"); return(NULL) }

  # Fixme: minor stochasticity here, perhaps due to numerical limitations?
  # Method indicates which test will be used; the higher the better score
  if (method == "hypergeometric") {
    enr <- enrichment.score(qofz, which.mode = response, annotation.samples = annotation.sample, method = method)
  }

  # This could be implemented, not sure how useful it would be
  # P(r|S) = P(S,r)/P(S) this assumes that all samples in S come from exactly one of the responses
  # prS <- P.rS(samples, model, pars = NULL, subnet.id, log = FALSE)
  # Or normalized version of the above: P(r|S)/P(r)
  
  if (method == "dependency") {
          
    # log(P(s,r)/P(s)P(r))
	  	    
    # now with actual sample density (not density mass as above)      
    # P(S) = sum_r P(S,r)

    #ps.log <- log(sum(get.P.rs.joint(s, models, subnet.id, log = FALSE)))

    # this requires features x samples matrix
    #psr.log <- P.s.r(t(dat), pars, log = TRUE)
    
    # log(P(s,r)/P(s)P(r)) = log(P(s|r)/P(s))
    #enr <- list(score = psr.log - ps.log, info = NULL)
   
  }
      
  if (method == "precision") {

    # This differs from 'hypergeometric' and 'dependency' in that they
    # compare proportion of factor level in response to factor level
    # in overall model now we compare proportion of factor level in
    # response to overall response (all samples in the response) This
    # does not necessarily correlate with the two other measures.  In
    # a way, this measures the purity of the response w.r.t. given
    # factor level
 
    # precision: TP/(TP + FP) = TP/n fraction of true posivites in
    # response -> here: density mass associated with this sample in
    # each response additionally normalize by the analogous fraction
    # in overall model density mass is the sum of individual sample
    # densities
             
    # P(s|r) / P(S|r)

    # density for each data point
    #dens <- sample.densities(s.ann, model, subnet.id, log = FALSE, summarize = FALSE)[response, s.ann]

    # P(s,r)/P(s)P(r) = P(s|r)/P(s) for factor level samples
    # Fraction of total density mass of factor level sample compared to all samples within the response
    # and w.r.t. overall density mass of the sample

    # relative density of sample
    #enr <- list(score = sum(dens[s])/sum(dens), info = NULL)

  }

  # recall: TP/(TP + FN) = TP/P fraction of all true positives
  # included in the response -> Add later

  # later utilize probabilistic interpretation of precision/recall? 
  
  enr

}






#' Enrichment for a specified sample group in the given response.
#' 
#' Calculate enrichment values for a specified sample group in the given
#' response.
#'
#' @param assignment.matrix Samples-to-modes assignment matrix. Soft or hard.
#' @param which.mode Mode for which the enrichment is calculated.
#' @param annotation.samples Samples assigned to the class for which the enrichment is calculated
#' @param method Enrichment method.
#' @param data Optionally required by some methods: data (samples x features)
#'
#' @return Enrichment score and info. Enrichment score varies from 0
#' (low) to 1 (high enrichment).
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse").
#' @keywords utilities
#' @export
#' @examples #

enrichment.score <- function (assignment.matrix, which.mode, annotation.samples, method = "hypergeometric", data = NULL) {

  # List all samples in the data set		 
  total.samples <- rownames(assignment.matrix)

  # List samples in the investigated mode
  subset.samples <- names(which(apply(assignment.matrix, 1, which.max) == which.mode))

  if (method == "hypergeometric") {

      N <- length(total.samples)

      # number of white balls in the urn
      m <- length(annotation.samples) 
    
      # number of black balls in the urn
      n <- N - m
      
      # number of balls drawn from the urn
      k <- length(subset.samples)  
      
      # overlap between investigated sample group among response
      # samples (using hard assignments)
      # number of white balls drawn without replacement
      q <- sum(annotation.samples %in% subset.samples)  
      
      # hypergeometric enrichment (small p, high enrichment)
      # take 1-p to indicate high enrichment with high score
      # use q-1 since lower.tail = FALSE indicates X > x calculation,
      # but we need X >=x
      # enr <- 1 - phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)
      pval <- phyper(q - 1, m, n, k, lower.tail = FALSE, log.p = FALSE)

      temp <- c(sample.size.total = N,
      	     	sample.size.subset = k, 
	     	sample.size.annotation = m,
	     	annotated.in.subset = q, 
		fraction.in.data = m/N,
		fraction.in.subset = q/k, 
		pvalue = pval)

      enr <- list(score = 1 - pval, info = temp)
  }

  enr

}


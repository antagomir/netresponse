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

#' Enrichment for a specified sample group in the given response.
#' 
#' Calculate enrichment values for a specified sample group in the given
#' response.
#'
#' @param subnet.id Subnet.
#' @param model NetResponseModel object.
#' @param s User-defined sample group. For instance, samples belonging to a particular annotation class.
#' @param response Response id (integer) within the subnet.
#' @param method Enrichment method.
#' @param data data (samples x features)
#'
#' @return List with enrichment statistics, depending on enrichment method.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso order.responses
#' @references See citation("netresponse").
#' @keywords utilities
#' @export
#' @examples #enr <- response.enrichment(subnet.id, model, sample, response, method)
#' 
response.enrichment <- function (subnet.id = NULL, model, s, response, method = "hypergeometric", data = NULL) {

  # s <- sample; method = "hypergeometric"; data = NULL

  if (is.null(data)) {
    data <- model@datamatrix
  }

  # samples x features
  if(is.vector(data)) {
    data2 <- matrix(data)
    rownames(data2) <- names(data)
    data <- data2
  }

  # pick sample data for the response and
  # ensure this is a matrix also when a single sample is given
  if (any(!s %in% rownames(data))) {
    warning("Not all samples are in the original data matrix; these are removed from enrichment analysis.")
    s <- intersect(s, rownames(data))
    s.ann <- s
  }

  if (class(model) == "NetResponseModel") {

    if (is.null(data)) {
      data <- model@datamatrix
    }

    if (is.numeric(subnet.id)) {
      subnet.id <- paste("Subnet", subnet.id, sep = "-")
      warning("subnet.id given as numeric; converting to character: ", subnet.id, sep="")
    }
 
     response.samples <- response2sample(model, subnet.id, component.list = TRUE)

     # Subnetwork feature names
     nodes <- model@subnets[[subnet.id]]
     dat <- matrix(data[, nodes, drop = FALSE], ncol = length(nodes))
     rownames(dat) <- rownames(data)
     colnames(dat) <- nodes
     # dat is now samples x features matrix

     pars <- get.model.parameters(model, subnet.id)
     dat <- model@datamatrix

   } else {

     # For mixture.model output
     response.samples <- response2sample(model, component.list = TRUE)
     pars <- model
     dat <- data
   }

  if (length(response.samples) == 0) { return(NULL) }
  response.sample <- response.samples[[response]]
  
  # Fixme: there is some minor stochasticity here, perhaps due to numerical limitations?
 
  # Method indicates which test will be used
  # FIXME: add other methods; the higher the better

  if (method == "hypergeometric") {

      N <- nrow(dat)

      # number of white balls in the urn
      m <- length(s) 
    
      # number of black balls in the urn
      n <- N - m
      
      # number of balls drawn from the urn
      k <- length(response.sample)  
      
      # overlap between investigated sample group among response samples (using hard assignments)
      q <- sum(s %in% response.sample)  #number of white balls drawn without replacement
      
      # hypergeometric enrichment (small p, high enrichment)
      # take 1-p to indicate high enrichment with high score
      # use q-1 since lower.tail = FALSE indicates X > x calculation, but we need X >=x
      # enr <- 1 - phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)
      pval <- phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)

      temp <- c(sample.size.total = N,
      	     			   sample.size.response = k, 
	     			   sample.size.mysample = m,
	     			   mysamples.in.response = q, 
				   fraction.in.data = m/N,
				   fraction.in.response = q/k, 
				   pvalue = pval)

      enr <- list(score = 1 - pval, info = temp)

  }
  
  # This could be implemented, not sure how useful it would be
  # P(r|S) = P(S,r)/P(S) this assumes that all samples in S come from exactly one of the responses
  # prS <- P.rS(samples, model, pars = NULL, subnet.id, log = FALSE)
  # Or normalized version of the above: P(r|S)/P(r)
  
  if (method == "dependency") {
          
    # log(P(s,r)/P(s)P(r))
	  	    
    # now with actual sample density (not density mass as above)      
    # P(S) = sum_r P(S,r)
    ps.log <- log(sum(get.P.rs.joint(s, model, subnet.id, log = FALSE)))

    # this requires features x samples matrix
    psr.log <- P.s.r(t(dat), pars, log = TRUE)
    
    # log(P(s,r)/P(s)P(r)) = log(P(s|r)/P(s))
    enr <- list(score = psr.log - ps.log, info = NULL)
   
  }
      
  if (method == "precision") {

    # This differs from 'hypergeometric' and 'dependency' in that they
    # compare proportion of factor level in response to factor level in overall model
    # now we compare proportion of factor level in response to overall response (all samples in the response)
    # This does not necessarily correlate with the two other measures.
    # In a way, this measures the purity of the response w.r.t. given factor level
 
    # precision: TP/(TP + FP) = TP/n fraction of true posivites in response
    # -> here: density mass associated with this sample in each response
    # additionally normalize by the analogous fraction in overall model
    # density mass is the sum of individual sample densities
             
    # P(s|r) / P(S|r)

    # density for each data point
    dens <- sample.densities(s.ann, model, subnet.id, log = FALSE, summarize = FALSE)[response, s.ann]

    # P(s,r)/P(s)P(r) = P(s|r)/P(s) for factor level samples
    # Fraction of total density mass of factor level sample compared to all samples within the response
    # and w.r.t. overall density mass of the sample

    # relative density of sample
    enr <- list(score = sum(dens[s])/sum(dens), info = NULL)

  }

  # recall: TP/(TP + FN) = TP/P fraction of all true positives included in the response
  # -> Not needed right now

  # later utilize probabilistic interpretation of precision/recall? 
  
  enr
}



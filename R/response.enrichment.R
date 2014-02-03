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
#' @param total.samples All samples in the data
#' @param response.samples Samples in the investigated subset
#' @param annotated.samples Samples at the investigated annotation level for enrichment calculation
#' @param method Enrichment method.
#'
#' @return List with enrichment statistics, depending on enrichment method.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso order.responses
#' @references See citation("netresponse")
#' @keywords utilities
#' @export
#' @examples #enr <- response.enrichment(subnet.id, models, sample, response, method)

response.enrichment <- function (total.samples, response.samples, annotated.samples, method = "hypergeometric") {

  if (length(response.samples) == 0) { warning("No samples in response"); return(NULL) }
  # Fixme: minor stochasticity here, perhaps due to numerical limitations?
  # Method indicates which test will be used; the higher the better score
  if (method == "hypergeometric") {
    enr <- enrichment.score(total.samples, response.samples, annotated.samples, method = method)
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







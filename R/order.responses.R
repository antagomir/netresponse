
# Copyright (C) 2010-2011 Leo Lahti
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

# Fixme: finish this later
#order.samples <- function (subnet.id, model, phenodata, which.factor, response, method = "hypergeometric") {
#    
#  # - for given response, order factor levels by association strength (enrichment score)
#  #   P(s|r) = P(s,r)/P(r) total sample density and/or average sample density (both for individuals and groups s/S)
#  #   P(s|r)/P(s) = P(s,r)/P(r)P(s) 
#
#  #   P(s|r)/P(s) -> OK: method = "dependency"
#  #   P(s|r) = P(s,r)/P(r)  -> OK, normalized version P(s|r)/P(S|r) available: method = "precision"
#
#  enr <- response.enrichments(subnet.id, model, phenodata, which.factor, response, method)
# 
#  # For the given response, return levels of the given factor (decreasing ordering by enrichement score)
#  sort(enr, decreasing = TRUE)
#}



order.responses <- function (model, sample, method = "hypergeometric", min.size = 2, max.size = Inf, min.responses = 2) {
  # Given sample (for instance set of samples associated with a given factor level)
  # order the responses across all subnetworks based on their association strength

  # For a given sample, calculate enrichment values in each response
  subnets <- responses <- scores <- pvals <- c()
  enrichment.info <- list()
  cnt <- 0

  # Get model statistics
  stat <- model.stats(model)

  # Filter the results
  sn <- get.subnets(model, get.names = TRUE, min.size, max.size, min.responses)
  stat <- stat[names(sn),]

  # Check enrichment in the selected responses  
  for (subnet.id in rownames(stat)) {

    for (response in 1:length(model@models[[subnet.id]]$w)) {
  
      enr <- response.enrichment(subnet.id, model, sample, response, method)

      # add further info about enrichments
      cnt <- cnt + 1
      enrichment.info[[cnt]] <- c(subnet = subnet.id, response = response, enrichment.score = enr$score, enr$info) 

    }
  }
			
  enr <- as.data.frame(t(sapply(enrichment.info, identity)))

  if ("pvalue" %in% colnames(enr)) {
    # calculate q-values
    library(qvalue)
    enr$qvalue <- qvalue(as.numeric(as.character(enr$pvalue)))$qvalues
  }

  enr[,3:ncol(enr)] <- apply(enr[,3:ncol(enr)], 2, as.numeric)
  ord <- order(enr$enrichment.score, decreasing = TRUE)
  enr <- enr[ord,]
  enr[["subnet"]] <- as.character(enr[["subnet"]])

  # Add subnet info in the result table
  enr <- cbind(enr, stat[enr$subnet,])

  list(ordered.responses = enr, method = method, sample = sample)
		
}
			    
response.enrichment <- function (subnet.id, model, s, response, method = "hypergeometric") {

  # s:   # samples associated with this factor level (ensure they are in the data)

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", subnet.id, sep="")
  }

  response.sample <- response2sample(model, subnet.id, component.list = TRUE)[[response]]

  # Fixme: there is some minor stochasticity here, perhaps due to numerical limitations?
 
  pars <- get.model.parameters(model, subnet.id)

  # All samples
  s.ann <- rownames(model@datamatrix) # model@samples

  # Subnetwork feature names
  nodes <- pars$nodes

  # pick sample data for the response and
  # ensure this is a matrix also when a single sample is given
  dat <- matrix(model@datamatrix[s, nodes], ncol = length(nodes))
  rownames(dat) <- s
  colnames(dat) <- nodes
  # dat is now samples x features matrix
      
  # Method indicates which test will be used
  # FIXME: add other methods; the higher the better

  if (method == "hypergeometric") {

      N <- nrow(model@datamatrix)

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
      #enr <- 1 - phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)
      pval <- phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)

      temp <- c(sample.size.total = N,
      	     			   sample.size.response = k, 
	     			   sample.size.mysample = m,
	     			   mysamples.in.response = q, 
				   fraction.in.data = m/N,
				   fraction.in.response = q/k, pvalue = pval)

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
    psr.log <- P.sr(t(dat), pars, log = TRUE)
    
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


list.significant.responses <- function (model, sample, qth = 1, method = "hypergeometric") {

  # Order responses according to their association with the given sample group
  o <- order.responses(model, sample = sample, method = "hypergeometric")$ordered.responses
  o[which(o$qvalue < qth),]

}



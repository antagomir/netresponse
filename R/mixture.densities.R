# Copyright (C) 2008-2012 Leo Lahti
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

###################################################################


 #  "An article about computational science in a scientific publication
 #  is not the scholarship itself, it is merely advertising of the
 #  scholarship. The actual scholarship is the complete software
 #  development environment and the complete set of instructions
 #  which generate the figures." 
 #                                          - Jon Claerbout

###################################################################

P.rS <- function (dat, pars, log = TRUE) {

  # Probability of a response, given sample (group)
  # P(r|S) = P(S|r)P(r)/P(S) = P(S, r)/(sum_r P(S, r))
  # NOTE: unstable due to overflows in particular when multiple samples are used
  
  #psr <- get.P.rs.joint(sample, model, subnet.id, log = FALSE)  
  # joint density P(S, r) for each component r
  psr <- P.rs.joint(dat, pars, log = FALSE)
  
  # Log P(r|S) = logP(S, r) - log(sumr(P(S, r)))
  logp <- log(psr) - log(sum(psr))
  
  if (log) {
    logp
  } else {
    exp(logp)
  }
  
}


P.r.s <- function (dat, pars, log = TRUE) {
  # P(r|s) for each response r and sample s

  qofz <- t(apply(P.rs.joint.individual(dat, pars, log = FALSE), 2, function (x) {x/sum(x)}))

  if ( log ) { qofz <- log(qofz) }

  matrix(qofz, ncol(dat))
}


P.rs.joint.individual <- function (dat, pars, log = TRUE) {

  # P(r,s) = P(s|r)P(r) = P(r|s)P(s)
  # P(r, s) for all samples and responses. Should hold sum_r P(r, s) = P(s) -> OK
  #  Alternatively: logp.joint <- t(t(prs.log) + ps.log)
  
  # FIXME: merge with P.rs.joint and/or P.rS to avoid redundancy
  
  psr.log <- P.sr(dat, pars, log = TRUE)
  pr.log <- log(pars$w)

  #Prs.log <- get.P.rs(model, subnet.id, log = TRUE)
  #ps <- sum_r P(s,r) = sum_r P(s|r)P(r) 
  # P(r,s) = P(s|r)P(r) = P(r|s)P(s)
  #ps <- colSums(exp(psr.log)*exp(pr.log)) # for each sample, sum over responses
  #ps.log <- log(ps)

  logp.joint <- psr.log + pr.log
  colnames(logp.joint) <- colnames(dat)
  
  if (log) {
    logp.joint
  } else {
    exp(logp.joint)
  }
  
}



P.rs.joint <- function (dat, pars, log = TRUE) {

  # Joint (Log) density of a given sample group: P(S, r) = P(S|r)P(r)
  # and each mixture component r.
  
  # First calculate
  # P(S|r) for each response; length of output equals to number of responses
  # i.e. logsum of the individual sample densities
  pSr.log <- P.Sr(dat, pars, log = TRUE) #get.P.Sr(sample, model, subnet.id, log = TRUE)  
  # Density for each mixture component
  pr.log <- log(pars$w)

  logp.joint <- pSr.log + pr.log

  # Alternative:
  # logp.joint <- rowSums(get.P.rs.joint.individual(sample, model, pars, subnet.id, log = TRUE))

  if (log) {
    logp.joint
  } else {
    exp(logp.joint)
  }  
  
}





P.Sr <- function (dat, pars, log = TRUE) {

  # dat is features x samples matrix
  psr <- P.sr(dat, pars, log = TRUE)

  # Returns responses x samples matrix
  # for each response, calculate logsum over samples
  if (log) {
     rowSums(psr)
  } else {
     exp(rowSums(psr))
  } 
}


P.sr <- function (dat, pars, log = TRUE) {
 
  # dat: features x samples matrix

  # Log probability density on each data point for each response
  # P(s|r)
    
  if (!nrow(dat) == ncol(pars$mu)) { stop("Dimensions in dat and pars do not match!") }

  # FIXME: in many cases density needs to be calculated just within a
  # single response, here calculated for all responses. Speedup by
  # having this for given response only.

  # subnetid <- "Subnet-2"
  # dat <- t(model@datamatrix[subnets[[subnetid]], ]) # all samples, given subnet
  # Psr(dat, pars)

  if ( is.vector(dat) ) { dat <- as.matrix(dat, nrow = length(dat)) }

  # responses x samples matrix P(s|r)
  psr <- matrix(NA, nrow = length(pars$w), ncol = ncol(dat)) 
  if (!is.null(colnames(dat))) { colnames(psr) <- colnames(dat) }
  
  for ( response in 1:length( pars$w ) ) {
    # Given the diagonal covariances, the density is product (log-sum)
    # over the densities for individual features (on each data point)
    psr[response, ] <- colSums(dnorm(dat, mean = as.numeric(pars$mu[response, ]), sd = as.numeric(pars$sd[response, ]), log = TRUE))		  
  }

  logp <- psr # responses x samples
  rownames(logp) <- names(pars$w)
  colnames(logp) <- colnames(dat)

  if (log) {
     logp
  } else {
    exp(logp)
  }

}


P.S <- function (dat, pars, log = TRUE) {

  # psi: individual sample densities P(s) given in _log domain_

  # FIXME: develop this s.t. calculates individual sample densities independently
  
  # Overall probability of sample s, given the Gaussian mixture model
  # P(s) = sum_r P(s, r)
  # ps <- log(sum(get.P.rs.joint(sample, model, subnet.id, log = FALSE)))

  # FIXME: numerically does not hold tightly that P(S) = sumr P(S, r) = sumr P(S|r)P(r)
  # the latter equality holds, P(S) is problematic. Differences are not big for examples
  # I checked, but they are still notable. Check in more detail this one.
  #sum(get.P.rs.joint(s, model, subnet.id, log = FALSE))
  #sum(get.P.Sr(s, model, pars = NULL, subnet.id, log = FALSE) * pars$w)

  ps <- sum(P.s.individual(dat, pars, log = TRUE)) # log sum
  
  if (log) {
    ps
  } else {
    exp(ps)
  }
}



P.s.individual <- function (dat, pars, log = TRUE) {

  # FIXME: merge with P.S (?)
  
  # Overall probability of sample s, given the model. 
  # individually for each sample

  # responses x samples
  # for each sample (column), density mass is the sum over joint densities on individual responses
  # P(s) = sum_r P(s, r) = sum_r P(s,r) = sum_r P(s|r)P(r)

  ps <- colSums(P.rs.joint.individual(dat, pars, log = FALSE))

  # two alternatives to calculate P(s)
  # log(colSums(get.P.rs.joint.individual(rsample, model, subnet.id, log = FALSE)))
  # log(sum(get.P.Sr(rsample, model, subnet.id, log = FALSE) * get.P.r(model, subnet.id, log = FALSE)))
  # get.P.s(rsample, model, subnet.id, log = TRUE) 

  # log only after summation:
  if (log) {
    log(ps)
  } else {
    ps
  }
}




###########################################################



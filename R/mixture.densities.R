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

#' Description: Probabiity of mode given multiple samples (ie. data matrix)
#' Mainly for internal use; documentation will be provided later. Tools for calculating densities with Gaussian mixture models.
#'
#' Arguments:
#'  @param dat features x samples data matrix for mixture modeling
#'  @param pars Gaussian mixture model parameters (diagonal covariances); list with elements mu (mean vectors), sd (covariance diagonals), w (weights). The mu and sd are component x features matrices, w is vector giving weight for each component.
#'  @param log Logical. Return densities in log domain.
#'
#' Returns:
#'   @return Probability density
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal utilities

P.rS <- function (dat, pars, log = TRUE) {

  # Probability of a response, given sample (group)
  # P(r|S) = P(S|r)P(r)/P(S) = P(S, r)/(sum_r P(S, r))
  # P(r|S) = P(S|r)P(r)/sum_r(P(S|r)P(r)) = (a*P(S|r)P(r))/(a*sum_r(P(S|r)P(r)))
  # -> log P(r|S) = log((a*P(S|r)P(r))/(a*sum_r(P(S|r)P(r)))) = log(a) + log(P(S|r)P(r)) - log(a) - log(sum_r(P(S|r)P(r)))
  # = log(a*P(S, r)) - log(sum_r(a*P(S, r)))
  # = log(a*psr) - log(sum_r(a*psr)) = log(a) + log(psr) - log(sum_r(exp(log(a) + log(psr))))
  # NOTE: unstable due to overflows in particular when multiple samples are used
  # Log P(r|S) = logP(S, r) - log(sumr(P(S, r)))
  #logp <- psr.log - log(sum(exp(psr.log)))
  # psr <- get.P.rs.joint(sample, model, subnet.id, log = FALSE)  

  # joint density P(S, r) for each component r
  psr.log <- P.rs.joint(dat, pars, log = TRUE)
  
  # density P(r|S), avoiding numerical overflows with log.a trick
  log.a <- 10 - max(psr.log)
  logp <- log.a + psr.log - log(sum(exp(log.a + psr.log)))

  # Compared output in univariate case to the direct calculation:
  #ps <- pars$w * dnorm(4.467782, mean = pars$mu, sd = pars$sd); ps/sum(ps)
  #P.rS(matrix(4.467782), pars, log = FALSE)
  # -> OK

  if (log) {
    logp
  } else {
    exp(logp)
  }
  
}



#' Description: Probabiity of mode given a sample (a data vector)
#' Mainly for internal use; documentation will be provided later. Tools for calculating densities with Gaussian mixture models.
#'
#' Arguments:
#'  @param dat features x samples data matrix for mixture modeling
#'  @param pars Gaussian mixture model parameters (diagonal covariances); list with elements mu (mean vectors), sd (covariance diagonals), w (weights). The mu and sd are component x features matrices, w is vector giving weight for each component.
#'  @param log Logical. Return densities in log domain.
#'
#' Returns:
#'   @return Probability density
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal utilities

P.r.s <- function (dat, pars, log = TRUE) {
  # P(r|s) for each response r and sample s

  qofz <- t(apply(P.rs.joint.individual(dat, pars, log = FALSE), 2, function (x) {x/sum(x)}))

  if ( log ) { qofz <- log(qofz) }

  matrix(qofz, ncol(dat))
}


#' Description: Joint probabiity density for mode and sample
#' Mainly for internal use; documentation will be provided later. Tools for calculating densities with Gaussian mixture models.
#'
#' Arguments:
#'  @param dat features x samples data matrix for mixture modeling
#'  @param pars Gaussian mixture model parameters (diagonal covariances); list with elements mu (mean vectors), sd (covariance diagonals), w (weights). The mu and sd are component x features matrices, w is vector giving weight for each component.
#'  @param log Logical. Return densities in log domain.
#'
#' Returns:
#'   @return Probability density
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal utilities

P.rs.joint.individual <- function (dat, pars, log = TRUE) {

  # P(r,s) = P(s|r)P(r) = P(r|s)P(s)
  # P(r, s) for all samples and responses. Should hold sum_r P(r, s) = P(s) -> OK
  #  Alternatively: logp.joint <- t(t(prs.log) + ps.log)
  
  # FIXME: merge with P.rs.joint and/or P.rS to avoid redundancy
  
  psr.log <- P.s.r(dat, pars, log = TRUE)
  pr.log <- as.vector(log(pars$w))

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



#' Description: Joint probabiity density for mode and sample group
#' Mainly for internal use; documentation will be provided later. Tools for calculating densities with Gaussian mixture models.
#'
#' Arguments:
#'  @param dat features x samples data matrix for mixture modeling
#'  @param pars Gaussian mixture model parameters (diagonal covariances); list with elements mu (mean vectors), sd (covariance diagonals), w (weights). The mu and sd are component x features matrices, w is vector giving weight for each component.
#'  @param log Logical. Return densities in log domain.
#'
#' Returns:
#'   @return Probability density
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal utilities

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
    return(logp.joint)
  } else {
    return(exp(logp.joint))
  }  
  
}





#' Description: Probabiity density for sample group given mode
#' Mainly for internal use; documentation will be provided later. Tools for calculating densities with Gaussian mixture models.
#'
#' Arguments:
#'  @param dat features x samples data matrix for mixture modeling
#'  @param pars Gaussian mixture model parameters (diagonal covariances); list with elements mu (mean vectors), sd (covariance diagonals), w (weights). The mu and sd are component x features matrices, w is vector giving weight for each component.
#'  @param log Logical. Return densities in log domain.
#'
#' Returns:
#'   @return Probability density
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal utilities

P.Sr <- function (dat, pars, log = TRUE) {

  # dat is features x samples matrix
  psr <- P.s.r(dat, pars, log = TRUE)

  # Returns responses x samples matrix
  # for each response, calculate logsum over samples
  if (log) {
     rowSums(psr)
  } else {
     exp(rowSums(psr))
  } 
}



#' Description: Probabiity density for sample given mode
#' Mainly for internal use; documentation will be provided later. Tools for calculating densities with Gaussian mixture models.
#'
#' Arguments:
#'  @param dat features x samples data matrix for mixture modeling
#'  @param pars Gaussian mixture model parameters (diagonal covariances); list with elements mu (mean vectors), sd (covariance diagonals), w (weights). The mu and sd are component x features matrices, w is vector giving weight for each component.
#'  @param log Logical. Return densities in log domain.
#'
#' Returns:
#'   @return Probability density
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal utilities

P.s.r <- function (dat, pars, log = TRUE) {

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




#' Description: Probabiity density for sample 
#' Mainly for internal use; documentation will be provided later. Tools for calculating densities with Gaussian mixture models.
#'
#' Arguments:
#'  @param dat features x samples data matrix for mixture modeling
#'  @param pars Gaussian mixture model parameters (diagonal covariances); list with elements mu (mean vectors), sd (covariance diagonals), w (weights). The mu and sd are component x features matrices, w is vector giving weight for each component.
#'  @param log Logical. Return densities in log domain.
#'
#' Returns:
#'   @return Probability density
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal utilities

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




#' Description: Probabiity density for individual sample 
#' Mainly for internal use; documentation will be provided later. Tools for calculating densities with Gaussian mixture models.
#'
#' Arguments:
#'  @param dat features x samples data matrix for mixture modeling
#'  @param pars Gaussian mixture model parameters (diagonal covariances); list with elements mu (mean vectors), sd (covariance diagonals), w (weights). The mu and sd are component x features matrices, w is vector giving weight for each component.
#'  @param log Logical. Return densities in log domain.
#'
#' Returns:
#'   @return Probability density
#'            
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords internal utilities

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





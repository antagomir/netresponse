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


info.criterion <- function (nparams, nlog, logp, criterion = "BIC") {

  # Calculate various information criteria

  if (criterion == "AIC") {         # Akaike IC 
    return(AIC(nparams, nlog, logp))
  } else if (criterion == "BIC") {  # Bayesian IC   
    return(BIC(nparams, nlog, logp))
  } else if (criterion == "AICc") { # Akaike for linear models, finite sample
    return(AICc(nparams, nlog, logp))
  }

}




BIC <- function (nparams, nlog, logp) {

  # Calculate Bayesian Information Criterion (BIC)

  # NOTE:
  # Original formulation (Schwartz) assumed that data is iid and likelihood is in exponential family.
  # However, BIC is later derived with considerably loosened assumptions, see e.g.
  # Cavanaugh et al.: Generalizing the Derivation of the Schwarz Information Criterion
  # it seems that the assumptions listed in Section 3 will quarantee the validity of BIC
  # for mixtures of exponential family distributions; confirm.
  
  # negative free energy is lower bound for log(P(D|H))
  # logp = -cost

  nparams*nlog - 2*logp
  
}

AIC <- function (nparams, nlog, logp) {

  # Calculate Akaike Information Criterion (AIC)
  # Note: nlog not used here but included for compatibility
  
  # negative free energy is lower bound for log(P(D|H))
  # logp = -cost

  2*(nparams - logp)
  
}

AICc <- function (nparams, nlog, logp) {

  # Calculate Akaike Information Criterion correction for finite
  # sample size (AICc)
  
  #AICc is AIC with a correction for finite sample sizes, giving a
  #greater penalty for extra parameters. Burnham & Anderson (2002)
  #strongly recommend using AICc, rather than AIC, if n is small or k
  #is large. Since AICc converges to AIC as n gets large, AICc
  #generally should be employed. Using AIC, instead of AICc, when n is
  #not many times larger than k2, increases the probability of
  #selecting models that have too many parameters, i.e. of overfitting.
  #The probability of AIC overfitting can be substantial, in some
  #cases. Brockwell & Davis (p. 273) advise using AICc as the primary
  #criterion in selecting the orders of an ARMA model for time series.
  #McQuarrie & Tsai ground their high opinion of AICc on extensive
  #simulation work with regression and time series. AICc was first
  #proposed by Hurvich & Tsai (1989). Different derivations of it are
  #given by Brockwell & Davis, Burnham & Anderson, and Cavanaugh. All
  #the derivations assume a univariate linear model with
  #normally-distributed errors (conditional upon regressors); if that
  #assumption does not hold, then the formula for AICc will usually
  #change. Further discussion of this, with examples of other
  #assumptions, is given by Burnham & Anderson (2002, ch.7). Note that
  #when all the models in the candidate set have the same k, then AICc
  #and AIC will give identical (relative) valuations. In that
  #situation, then, AIC can always be used.
  
  # negative free energy is lower bound for log(P(D|H))
  # logp = -cost

  n <- exp(nlog)
  
  AIC(nparams, nlog, logp) + 2*nparams(nparams + 1)/(n - nparams - 1)
  
}




#AIC.c <- cmpfun(AIC)
#AICc.c <- cmpfun(AICc)
#BIC.c <- cmpfun(BIC)

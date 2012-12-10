# Copyright (C) 2008-2012 Olli-Pekka Huovilainen and Leo Lahti 
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
#

# The best academic advice I ever got was: "Spend at least an hour
# every day on the manuscript closest to publication".

#' generate.toydata
#' 
#' @usage D <- generate.toydata()
#' @param Dim Dimensionality of data
#' @param Nc Number of modes
#' @param Ns Number of data points
#' @param sd0 Component spread
#' @param rgam.shape Shape parameter for Gamma distribution 
#' @param rgam.scale Scale parameter for Gamma distribution
#' @return Simulated data matrix (samples x features)
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @examples # D <- generate.toydata()

generate.toydata <- function (Dim = 3, Nc = 3, Ns = 200, sd0 = 3, rgam.shape = 2, rgam.scale = 2, rseed = 12346) {

  set.seed(rseed)

  # Generate means and variances (covariance diagonals) for the components 
  component.means <- matrix(rnorm(Nc*Dim, mean = 0, sd = sd0), nrow = Nc, ncol = Dim)
  component.vars <- matrix(1/rgamma(Nc*Dim, shape = rgam.shape, scale = rgam.scale), 
	                 nrow = Nc, ncol = Dim)
  component.sds <- sqrt(component.vars)

  # Size for each component -> sample randomly for each data point from uniform distr.
  # i.e. cluster assignments
  sample2comp <- sample.int(Nc, Ns, replace = TRUE)
  D <- array(NA, dim = c(Ns, Dim))
  for (i in 1:Ns)  {
    # component identity of this sample
    ci <- sample2comp[[i]]
    cm <- component.means[ci,]
    csd <- component.sds[ci,]
    D[i,] <- rnorm(Dim, mean = cm, sd = csd)
  }

  colnames(D) <- paste("Feature-",  1:ncol(D), sep = "")
  rownames(D) <- paste("Sample-", 1:nrow(D), sep = "")

  list(data = D, means = component.means, sds = component.sds, sample2comp = sample2comp)
}

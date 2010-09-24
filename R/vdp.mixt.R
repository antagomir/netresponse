  #
  #  This file is a part of the NetResponse R package.
  #
  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti.
  #  Contact: Leo Lahti <leo.lahti@iki.fi>
  #
  #  This program is free software; you can redistribute it and/or
  #  modify it under the terms of the GNU General Public License as
  #  published by the Free Software Foundation; either version 2 of
  #  the License, or (at your option) any later version.
  #
  #  This program is distributed in the hope that it will be useful,
  #  but WITHOUT ANY WARRANTY; without even the implied warranty of
  #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  #  General Public License for more details
  #  <http://www.gnu.org/licenses/>.
  #
  #  This file is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner
  #


vdp.mixt <-
function(dat,
    prior.alpha    = 1,  
    prior.alphaKsi = 0.01,
    prior.betaKsi  = 0.01,
           do.sort = TRUE,   # qOFz sorted in decreasing fashion based on colSums(qOFz) 
         threshold = 1.0e-5, # minimal free energy improvement that stops the algorithm
         initial.K = 1,      # initial number of components
               ite = Inf,    # used on updatePosterior: maximum number of iterations
         implicit.noise = 0, # Adds implicit noise in vdp.mk.log.lambda.so and vdp.mk.hp.posterior.so
             c.max = 10      # max. candidates to consider in find.best.splitting.
                     )
{

  
  # Prior parameters
  opts <- list(
    prior.alpha    = prior.alpha,
    prior.alphaKsi = prior.alphaKsi,
    prior.betaKsi  = prior.betaKsi,
           do.sort = do.sort,
         threshold = threshold,
         initial.K = initial.K,      
               ite = ite,    
  implicitnoisevar = implicit.noise,
             c.max = c.max
  )

                         data  <- list()
  data[["given.data"]]         <- list()
  data[["given.data"]][["X1"]] <- dat

  qOFz         <- rand.qOFz(nrow(dat), initial.K)
  hp.prior     <- mk.hp.prior(data, opts)
  hp.posterior <- mk.hp.posterior(data, qOFz, hp.prior, opts)

  # Note: greedy gives components in decreasing order by size
  templist     <- greedy(data, hp.posterior, hp.prior, opts)
  templist$hp.prior <- c(templist$hp.prior, list(qOFz = qOFz))
  qOFz <- matrix(templist$hp.posterior$qOFz, nrow(dat))

  # number of mixture components (nonempty components only!)
  Kreal <- sum(colSums(qOFz) > 0)
  qOFz  <- matrix(qOFz[, 1:Kreal], nrow(dat))
    
  ###############################################
  
  # Calculate mixture model parameters
  # FIXME: move this outside from this vdp.mixt function
  
  # Negative free energy is (variational) lower bound for P(D|H) Use
  # it to approximate P(D|HClist <- list(C))

  # number of parameters
  # (d-dim. centroid + diag. cov. matrix + component weight for each Kreal)
  Nparams <- Kreal*(2*ncol(dat)  + 1)

  # Calculate map estimates of model parameters from the posterior

  # variances are assumed inverse Gamma distributed and here beta/alpha gives the expectation)
  variances  <-  templist$hp.posterior$KsiBeta/templist$hp.posterior$KsiAlpha  # Cluster variances
  variances  <-  matrix(variances[1:Kreal,], Kreal)

  # Ignore empty components assuming that the components have been
  # ordered in decreasing order by size
  means     <- templist$hp.posterior$Mubar    # component centroids
  means    <-  matrix(means[1:Kreal,], Kreal)

  #############################################
  
  # FIXME: improve later, or retrieve weight directly from vdp code
  # estimate weights using data and other parameters
  # take average to get more robust estimate over data points

  ws  <- array(NA, dim = c(nrow(dat), Kreal))
  for (datapoint in 1:nrow(dat)) {
    ws[datapoint,] <- compute.weight(qOFz, means, variances, dat, datapoint)             
  }

  # Remove those with NAs (zero densities in some components mess the weights)
  w <- apply(ws, 2, function (wcol) { mean(na.omit(wcol)) })

 
  #######################################
  
  results <- list(
                  hp.prior     = templist$hp.prior,
                  opts         = opts,
                  free.energy  = templist$free.energy,
                  hp.posterior = templist$hp.posterior,
                  K = Kreal,
                  qOFz = qOFz,
                  Nparams = Nparams,
                  means = means,
                  variances = variances,
                  weights = w 
                  )
  results
}


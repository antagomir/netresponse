vdp.mixt <-
function(dat,
    prior.alpha    = 1,  
    prior.alphaKsi = 0.01,
    prior.betaKsi  = 0.01,
           do.sort = TRUE,   # qOFz sorted in decreasing fashion based on colSums(qOFz) 
         threshold = 1.0e-5, # minimal free energy improvement that stops the algorithm
         initial.K = 1,      # initial number of components
               ite = Inf,    # used on updatePosterior: maximum number of iterations
         implicit.noise = 0, # Adds implicit noise in 
	 		     # vdp.mk.log.lambda.so and vdp.mk.hp.posterior.so
             c.max = 10,     # max. candidates to consider in find.best.splitting. 
	                     # i.e. truncation parameter
                             # Candidates are chosen based on their Nc value 
                             # (larger = better). Nc = colSums(qOFz)
           speedup = TRUE,        
	   	             # speedup: during DP, components are splitted
			     # based on their first PCA component.
                             # To speed up, approximate by using only subset 
			     # data to calculate PCA.
         min.size = 5 # min size for a component to be splitted
         ) {

 
  #
  #  This file is a part of the NetResponse R package.
  #
  #  Copyright (C) 2008-2011 Antonio Gusmao and Leo Lahti.
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

# INPUT: mydat::
#       Each Row is an observation.
#       Each Column is a variable.
#
# OUTPUT:
#    list(hp_prior, opts, free_energy, hp_posterior, K);
#
#    * hp_prior: prior info
#         - hp_prior$q_of_z: prior on observation labels
#         - Mu_mu: centroids, 
#         - S2_mu: variance,
#
#    * opts: option list used in training
#
#    * free_energy: free energy of mixture model found.
#
#    * hp_posterior = templist$hp_posterior
#
#    * K: Number of mixture components (clusters)
#  
#
################  ALGORITHM SUMMARY  ################
# This code implements Gaussian mixture models with diagonal covariance matrices. 
# The following greedy iterative approach is taken in order to obtain the number
# of mixture models and their corresponding parameters:
#
# 1. Start from one cluster, $T = 1$.
# 2. Select a number of candidate clusters according to their values of 
#    "Nc" = \sum_{n=1}^N q_{z_n} (z_n = c) (larger is better).
# 3. For each of the candidate clusters, c: 
#     3a. Split c into two clusters, c1 and c2, through the bisector of its 
#         principal component. Initialise the responsibilities 
#         q_{z_n}(z_n = c_1) and q_{z_n}(z_n = c_2). 
#     3b. Update only the parameters of c1 and c2 using the observations that
#         belonged to c, and determine the new value for the free energy, F{T+1}.
#     3c. Reassign cluster labels so that cluster 1 corresponds to the largest 
#         cluster, cluster 2 to the second largest, and so on.
# 4. Select the split that lead to the maximal reduction of free energy, F{T+1}.
# 5. Update the posterior using the newly split data.
# 6. If FT - F{T+1} < \epsilon then halt, else set T := T +1 and go to step 2.
#
# The loop is implemented in the function greedy(...)

#   prior.alpha = 1; prior.alphaKsi = 0.01; prior.betaKsi = 0.01;
#    do.sort = TRUE; threshold = 1.0e-5; initial.K = 1; ite = Inf;
#    implicit.noise = 0; c.max = 10


  #system("/Rpath/bin/R CMD SHLIB /path/netresponse.c")
  #dyn.load("/path/netresponse/src/netresponse.so")

  # Prior parameters
  opts <- list(
    prior.alpha    = prior.alpha,    # Remark: result is quite insensitive to this variable
    prior.alphaKsi = prior.alphaKsi, # smaller -> less clusters (and big!) -> quite sensitive
    prior.betaKsi  = prior.betaKsi,  # larger -> less clusters (and big!)
           speedup = speedup,        # speed up calculations
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

  # sample-component assignments
  qOFz         <- rand.qOFz(nrow(dat), initial.K)

  # The hyperparameters of priors
  hp.prior     <- mk.hp.prior(data, opts)

  # Posterior
  hp.posterior <- mk.hp.posterior(data, qOFz, hp.prior, opts)

  # Note: greedy gives components in decreasing order by size
  templist     <- greedy(data, hp.posterior, hp.prior, opts, min.size)
  templist$hp.prior <- c(templist$hp.prior, list(qOFz = qOFz))
  qOFz <- matrix(templist$hp.posterior$qOFz, nrow(dat))

  ###############################################

  # Retrieve model parameters 

  # number of mixture components (nonempty components only!)
  # response must have 'non-negligible' probability mass!
  # i.e. at least some points associated with it
  Kreal <- max(apply(qOFz, 1, which.max))  #sum(colSums(qOFz) > 1e-3)
  qOFz  <- matrix(qOFz[, 1:Kreal], nrow(dat))
      
  # Calculate mixture model parameters
  # FIXME: move this outside from this vdp.mixt function
  
  # Negative free energy is (variational) lower bound for P(D|H) Use
  # it to approximate P(D|HClist <- list(C))

  # number of parameters
  # (d-dim. centroid + diag. cov. matrix + component weight for each Kreal)
  Nparams <- Kreal*(2*ncol(dat)  + 1)

  # Calculate map estimates of model parameters from the posterior

  # variances are assumed inverse Gamma distributed and here beta/alpha gives the expectation)

  # Parameters of the inverse Gamma function for component variances
  invgam.shape <- matrix(templist$hp.posterior$KsiAlpha[1:Kreal,], Kreal)
  invgam.scale <- matrix(templist$hp.posterior$KsiBeta[1:Kreal,],  Kreal)

  # Calculate variances (mean and mode of the invgam distr.) from scale and shape
  # FIXME: beta/alpha used in C code
  #var.update <- matrix(invgam.scale/invgam.shape, Kreal)
  #var.mean <- matrix(invgam.scale/(invgam.shape - 1), Kreal)
  var.mode  <- matrix(invgam.scale/(invgam.shape + 1), Kreal)
  variances <- var.mode # select mean, mode, or their average for update

  # Ignore empty components assuming that the components have been
  # ordered in decreasing order by size
  # component centroids
  centroids <- matrix(templist$hp.posterior$Mubar[1:Kreal,], Kreal)

  #############################################
  
  # FIXME: improve later, or retrieve weight directly from vdp code
  # estimate weights using data and other parameters
  # take average to get more robust estimate over data points

  # FIXME: can be sped up with apply. Need to modify compute.weight also
  #ws  <- array(NA, dim = c(nrow(dat), Kreal))
  #for (datapoint in 1:nrow(dat)) {
  #  ws[datapoint,] <- compute.weight(qOFz, centroids, variances, dat, datapoint)             
  #}
  # Calculate weights on at most 20 points
  # (weights are almost identical in all points, calculate on many points to increase robustness)
  rsample <- sample(nrow(dat), min(nrow(dat), 20))
  ws <- matrix(sapply(rsample, function (i) {compute.weight(qOFz[i,], centroids, variances, dat[i,])}), length(rsample))

  # Remove those with NAs (zero densities in some components mess the weights)
  w <- apply(ws, 2, function (wcol) { mean(wcol, rm.na = TRUE) })

  posterior <- list(
                  weights = w, 
                  centroids = centroids, # Mubar
                  sds = sqrt(variances), # alpha, beta
                  qOFz = qOFz,
		  Nc = colSums(qOFz), # component sizes
		  invgam.shape = invgam.shape, # KsiAlpha
		  invgam.scale = invgam.scale,  # KsiBeta
                  Nparams = Nparams, # number of model parameters
                  K = Kreal # number of components
  )

  # Later: include these from hp.posterior to the output later if needed.
  # "Mutilde[1:Kreal,]"  "gamma[,1:Kreal]"  "Uhat"

  list(prior = templist$hp.prior, 
       posterior    = posterior,  
       opts         = opts,  
       free.energy  = templist$free.energy)
}




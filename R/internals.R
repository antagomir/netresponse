get.subnet <- function (res, subnet.id) {

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", subnet.id, sep="")
  }
    
  # Nodes for a given subnet
  get.subnets(res)[[subnet.id]]

}


compute.weight <- function (pt, mu, vars, xt) {

  # pt <- qOFz[t,] # P(c|t) for all c
  # xt  <- dat[t, ] # data point
  
  # Initial weights are zero
  w <- rep.int(0, nrow(mu))  
  sds <- sqrt(vars)

  # Avoid overflows: no probability can be 0 exactly. Add negligible constant in these cases
  pt[pt < 1e-320] <- 1e-320
  pt <- pt/sum(pt) # renormalize

  # set arbitrary weigh for the first cluster
  # (only relations between weights matter)
  # normalize later
  #w[[1]] <- 1 # added below
  # operate on log domain to avoid floating errors

  logdens <- sum(dnorm(xt, mu[1, ], sds[1, ], log = TRUE))
 
  if (nrow(mu) > 1) {
    logw <- sapply(2:nrow(mu), function (i) { 
      #print(prod(dnorm(xt, mu[i, ], sds[i, ])))
  
      #dens * pt[[i]]/(pt[[1]] * prod(dnorm(xt, mu[i, ], sds[i, ]))) 

      logdens + log(pt[[i]]) - log(pt[[1]]) - sum(dnorm(xt, mu[i, ], sds[i, ], log = TRUE))

    })
    w <- exp(logw)

    # Densities in original domain, standardized s.t. first weight is 1
    w <- c(1, w)

    # normalize to unity and return
    return(w/sum(w))
  } else {
    return(1) #w <- 1
  }
    
}



############################################


retrieve.model <- function (model, subnet.id) {

  # Recalculate the model for a given subnet
  # Note: the algorithm has some stochasticity in initialization etc.
  # so the results may not be exactly same each time
  
  # Former: get.model

  # model: output from run.netresponse function
  # subnet.id: id/index of the subnet to check
  # level: which agglomeration step
  # datamatrix for which the model was calculated
  
  #  Copyright (C) 2008-2012 Leo Lahti
  #  Licence: GPL >=2

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", subnet.id, sep="")
  }
  
  # Get subnet nodes
  nodes <- model@subnets[[subnet.id]]

  # Compute the model
  x <- matrix(model@datamatrix[, nodes], nrow(model@datamatrix))

  pars <- mixture.model(x, mixture.method = model$params$mixture.method, 
  		   max.responses = model$params$max.responses,
		   implicit.noise = model$params$implicit.noise, 
		   prior.alpha = model$params$prior.alpha,  
		   prior.alphaKsi = model$params$prior.alphaKsi,
		   prior.betaKsi = model$params$prior.betaKsi,
		   vdp.threshold = model$params$vdp.threshold, 
		   initial.responses = model$params$initial.responses, 
		   ite = model$params$ite, 
		   speedup = model$params$speedup, 
		   bic.threshold = model$params$bic.threshold, 
		   pca.basis = model$params$pca.basis)


  pars
}



# INPUT:   data, hp_posterior, hp_prior, opts
# OUTPUT:  list(free_energy,hp_posterior,data,c)
#        * free_energy: free energy of the best split found
#        * hp_posterior: posterior info of the best split found
#        * c: index of the cluster that resulted in the best split found
#
# DESCRIPTION: Implements the VDP algorithm steps 2 to 4.

find.best.splitting <- function(data, hp.posterior, hp.prior, opts, min.size = 5){

  # min.size: minimum size of a component required for splitting
  
  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner  

  dat <- data$given.data$X1
  epsilon <- 1e-10    # FIXME: should this be a tunable function parameter?
  c.max <- opts$c.max 

  # ALGORITHM STEP 2
 
  # Sort clusters by size, use at most c.max candidates and ensure cluster size is > 2
  candidates <- order(hp.posterior$Nc, decreasing = TRUE)
  candidates <- candidates[hp.posterior$Nc[candidates] > 2]
  # candidates <- which(hp.posterior$Nc > 2)
  if ( length(candidates) == 0 ) { c <- 1 }

  qOFz <- mk.qOFz(data, hp.posterior, hp.prior, opts)
  fc   <- mk.E.log.q.p.eta(data, hp.posterior, hp.prior, opts)
  log.lambda <- mk.log.lambda(data, hp.posterior, hp.prior, opts) 
  sub.data  <- data

  # ALGORITHM STEP 3 (3a,3b,3c) check which split gives best improvements in cost

  new.qOFz.list <- list()   #Initialize
  new.free.energy <- rep(Inf, min(c.max, length(candidates)))

  for (c in candidates[1:min(c.max, length(candidates))]) {

    # relating.n has the indexes of data points that belonged to the candidate cluster prior
    # to splitting (that's why it is the sum over the now 2 clusters (after splitting).
    # REMARK: Is 0.5 ok? - when there are lots of clusters it is natural to 
    # assume some points will have less than 0.5 for any cluster.

    relating.n <- which(qOFz[, c] > 0.5)
    if (length(relating.n) == 0) { next } else {}

    # ALGORITHM STEP 3a. split the candidate cluster
    # Split cluster c in qOFz into two smaller ones

    new.c     <- ncol(qOFz)
    new.qOFz  <- split.qofz(qOFz, c, new.c, dat, opts$speedup, min.size)

    new.K     <- ncol( new.qOFz )
    sub.qOFz  <- new.qOFz[relating.n, unique(c(c, new.c, new.K))]
    
    # Ensure it remains a matrix
    sub.data$given.data$data <- sub.data$given.data$X1 <- array(dat[relating.n, ],
                                    dim = c(length(relating.n), ncol(dat)))    

    # ALGORITHM STEP 3b and 3c
    # update the posterior of the split clusters for a small number of iter.
    # update_posterior sorts clusters by size

    sub.hp.posterior <- mk.hp.posterior(sub.data, sub.qOFz, hp.prior, opts)
    dummylist        <- updatePosterior(sub.data, sub.hp.posterior, hp.prior, opts, 10, 0)
    sub.hp.posterior <- dummylist$hp.posterior
    sub.qOFz         <- dummylist$qOFz

    # FIXME: check this already previously for c == new.c? 
    if(ncol( sub.qOFz ) < 3) { next } else { } 

    # If there are more than 1 empty components then go to next step
    if(sum(colSums(sub.qOFz) < epsilon) > 1) { next } else { }

    sub.log.lambda <- mk.log.lambda(data, sub.hp.posterior, hp.prior, opts)
    insert.indices <- c(c, new.c, new.K:(new.K + ncol(sub.qOFz) - 3))
    
    if(max(insert.indices) > ncol(log.lambda)){
      new.log.lambda <- cbind(log.lambda, array(0, dim = c(nrow(log.lambda), max(insert.indices) - ncol(log.lambda))))
    } else{ new.log.lambda <- log.lambda }

    
    new.log.lambda[, insert.indices] <- sub.log.lambda
    new.fc <- fc
    new.fc[insert.indices] <- mk.E.log.q.p.eta(sub.data, sub.hp.posterior, hp.prior, opts)
    new.free.energy[[c]] <- mk.free.energy(data, sub.hp.posterior, hp.prior, opts, new.fc, new.log.lambda)$free.energy

    # if new.qOFz is not large enough to accommodate update from sub.qOFz then add columns
    if (ncol(new.qOFz) < max(insert.indices)) {
      new.qOFz <- cbind(new.qOFz, array(0, dim = c(nrow(new.qOFz), max(insert.indices) - ncol(new.qOFz))))
    }
    
    new.qOFz[relating.n, ] <- 0
    new.qOFz[relating.n, insert.indices] <- sub.qOFz
    new.qOFz.list[[c]]     <- new.qOFz
  }

  # Select cluster split that minimizes free energy
  c <- which.min(new.free.energy)
  free.energy <- new.free.energy[[c]]  

  
  if(is.infinite(free.energy)){
    c <- -1
  } else { hp.posterior <- mk.hp.posterior(data, new.qOFz.list[[c]], hp.prior, opts) }

  list( free.energy = free.energy,
       hp.posterior = hp.posterior,
               data = data,
                  c = c)
}




find.best.neighbor <- function (G, max.subnet.size, network, delta) {

  # Order edges by delta values. The two subnets with the smallest delta are
  # joined unless the merged subnet exceeds max size.
  o <- order(delta)

  # Check size of the resulting merged subnet one-by-one, starting from the smallest
  best.found <- FALSE
  cnt <- 0
  best.edge <- a <- b <- NULL
  mindelta <- Inf

  while (!best.found) {
    cnt <- cnt + 1
    ind <- o[[cnt]] # FIXME: could use the filtering for top neighborghs here as well
    z <- network[1, ind]
    i <- network[2, ind]    

    # Finish when the new merged subnetwork, which does not exceed max size, is found
    if (length(c(G[[z]], G[[i]])) <= max.subnet.size){
      # Sort a and b
      a <- min(c(z, i))
      b <- max(c(z, i))
      best.edge <- ind
      mindelta <- delta[[best.edge]]
      best.found <- TRUE
    } else {
      # Infinite cost for merges that exceed maximum size
      delta[[ind]] <- Inf 
    }
  }

  list(a = a, b = b, mindelta = mindelta, best.edge = best.edge, delta = delta)

}


join.subnets <- function (network, delta, best.edge) {

  # Pick the nodes to merge
  a <- network[1, best.edge]
  b <- network[2, best.edge]  
  
  # replace b nodes by a: this 
  # replaces edges from other nodes pointing to b 
  # to point into a
  network[network == b] <- a
  
  # remove delta values for nodes associated with node a
  # since this node now has different (merged) set of features
  # and the model needs to be recalculated
  inds <- apply(network == a, 2, any)
  delta[inds] <- NA

  # sort entries (row1 < row2)
  network <- matrix(apply(network, 2, sort), 2)
      
  list(network = network, delta = delta)

}






# INPUT: matrix (q_of_z)
# OUTPUT: matrix

# DESCRIPTION: Sorted matrix in decreasing fashion based on the value
#              of colSums.  Remark: The last column of the matrix is
#              kept in place (it is not sorted).

sortqofz <- function(qOFz){

  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner
  
  Nc <- colSums(qOFz)
  I  <- order(Nc[-length(Nc)], decreasing = TRUE) 
  I  <- c(I, ncol(qOFz)) # add last column element separately (outside of ordering)
                                    
  # Order the cols and ensure qOFz remains a matrix
  array(qOFz[, I], dim = dim(qOFz))
  
}


# INPUT:   data   - matrix with data vectors
#          K      - number of clusters
# OUTPUT:  q_of_z - matrix of size N*(K+1).
# DESCRIPTION: This function assigns data randomly to K clusters by drawing cluster
#              membership values from a uniform distribution.

rand.qOFz <- function(N, K){

  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  qOFz <- matrix( runif(N*(K+1)), N, K + 1)
  qOFz[, K + 1] <- 0

  # normalize and return (each row should sum up to 1)
  qOFz/rowSums(qOFz) 

}




# INPUT: "old" free_energy value, "new" free_energy value, options
# OUTPUT: bool: 0 if the there was no significant improvement.
#               1 if new_free_energy is smaller than free_energy (more than opts$threshold).



free.energy.improved <- function(free.energy, new.free.energy,
                                 warn.when.increasing, threshold)
{

  #INPUT: "old" free.energy value, "new" free.energy value, options
  #OUTPUT: bool: 0 if the there was no significant improvement.
  #              1 if new.free.energy is smaller than free.energy
  #                (more than opts$threshold).

  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  diff <- new.free.energy - free.energy

  v <- abs(diff/free.energy)
  
  if(is.nan(v) || v < threshold){
    bool <- 0
  } else { 
    if(diff > 0){
      if( warn.when.increasing ){
        if( v > 1e-3 ){
          stop(c("the free energy increased. The diff is ", toString(diff)))
        } else {
          warning(c("the free energy increased. The diff is ", toString(diff)))
        }
      }
      bool <- 0
    } else {
      bool <- ifelse(diff == 0, 0, 1)
    }
  }
  
  bool
}




# INPUT:   data: structure with data matrix
#          hp_posterior: posterior information for the current mixture model
#          hp_prior: prior information
#          opts: options list.
#
# OUTPUT:  list(free_energy, hp_posterior, hp_prior, data);
# 
# DESCRIPTION: Read the main description on the beginning of the file.


greedy <- function(data, hp.posterior, hp.prior, opts, min.size){

  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner
  
  free.energy <- mk.free.energy(data, hp.posterior, hp.prior, opts)$free.energy

  while(1){

    # ALGORITHM STEP 2-4

    templist <- find.best.splitting(data, hp.posterior, hp.prior, opts, min.size)
    new.free.energy  <- templist$free.energy
    new.hp.posterior <- templist$hp.posterior
    c                <- templist$c
    if ( c == (-1) ) { break } # infinite free energy -> break splitting

    # ALGORITHM STEP 5
    dummylist <- updatePosterior(data, new.hp.posterior, hp.prior,
                                  opts, ite = opts$ite, do.sort = 1)

    new.free.energy  <- dummylist$free.energy
    new.hp.posterior <- dummylist$hp.posterior
    
    # ALGORITHM STEP 6
    if( free.energy.improved(free.energy, new.free.energy, 0, opts$threshold) == 0 ) {
      break #free.energy didn't improve, greedy search is over
    } 

    free.energy  <- new.free.energy
    hp.posterior <- new.hp.posterior

  }

  list(free.energy = free.energy,
       hp.posterior = hp.posterior,
       hp.prior = hp.prior,
       data = data)
}

###############################################################################


#INPUT:   data, hp_posterior, hp_prior, opts, ite, do_sort
#             * do_sort: TRUE/FALSE: indicates whether the clusters should be 
#                        sorted by size or not.
#             * ite: number of update iterations.
#OUTPUT:  Updated parameters: list(free_energy, hp_posterior, q_of_z)
#DESCRIPTION: Updates the posterior of the mixture model. if do_sort=true it also
#             sorts the cluster labels by size. (i.e. cluster 1 = largest cluster)


updatePosterior <- function(data, hp.posterior, hp.prior, opts, ite = Inf, do.sort = 1) {
                            
  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  epsilon <- 1e-10
  free.energy <- Inf
  i <- last.Nc <- start.sort <- 0
  opts.internal <- opts
  break.loop <- FALSE
  
  while( 1 ){ 
    i <- i + 1

    new.free.energy <- Inf
    cnt <- 0
    while (is.infinite(new.free.energy) && cnt <= 10) {
      
      templist <- mk.free.energy(data, hp.posterior, hp.prior, opts.internal)
      new.free.energy <- templist$free.energy
      log.lambda <- templist$log.lambda
      if( is.infinite(new.free.energy) ) {
        warning("Free energy not finite: adding implicit noise.")
        opts.internal$implicitnoisevar <- opts.internal$implicitnoisevar + 0.1
      }
      cnt <- cnt + 1
    }

    if ( (is.finite(ite) && i>= ite) ||
         (is.infinite(ite) && 
          free.energy.improved(free.energy, new.free.energy, 0, opts.internal$threshold) == 0)){
        free.energy <- new.free.energy
        if( do.sort && opts.internal$do.sort && (!start.sort) && is.finite(free.energy)){
          start.sort <- 1
        } else break # this will break the while loop
      }

      last.Nc <- hp.posterior$Nc
      free.energy <- new.free.energy
      qOFz <- mk.qOFz(data, hp.posterior, hp.prior, opts.internal, log.lambda)

      # if the last component is not 'empty' and max number components not reached
      # add a new empty component
      if(sum(qOFz[, ncol(qOFz)]) >= epsilon && ncol(qOFz) < opts$c.max){ 
        qOFz <- cbind(qOFz, 0) 
      }

      # Sort components by size (note: last component kept in its place)
      if( start.sort ){ qOFz <- sortqofz(qOFz) }

      # Pick at most c.max+1 components, remove the smallest one, except the one added in this iteration (ie. c.max + 1)
      # FIXME: consider how to improve the implementation regarding to this!
      # (3/2012)
      #if (ncol(qOFz) == opts$c.max + 1) {
      # 	 qOFz <- qOFz[, setdiff(1:ncol(qOFz), opts$c.max), drop = FALSE]        
      #} 

      # If the smallest of the previous components
      # (excluding the one added in this interation)
      # is empty then remove it
      # (if new component was not added this is ok to have as well)
      if(sum(qOFz[, ncol(qOFz) - 1 ]) < epsilon){
        qOFz <- matrix(qOFz[, -(ncol(qOFz) - 1)], nrow(qOFz))
      }

      qOFz <- matrix(qOFz/rowSums(qOFz), nrow(qOFz)) # probabilities sum to one      
      hp.posterior <- mk.hp.posterior(data, qOFz, hp.prior, opts.internal)    

  }

  list( free.energy = free.energy,
       hp.posterior = hp.posterior,
               qOFz = qOFz)
}


sumlogsumexp <- function(log.lambda){.Call("vdpSumlogsumexp", log.lambda, PACKAGE = "netresponse")}


softmax <- function( A ){
  qOFz <- .Call("vdpSoftmax", A, PACKAGE = "netresponse")

  # Ensure that this remains a matrix even if it has 1-dimension on rows or cols
  qOFz <- array(qOFz, dim = dim(A))
  
  # In rare (< 1 / 1e6?) cases, and particularly when sample size is
  # large (>1500) the vdpSoftmax function seems to produce NaNs. This
  # occurred for particularly low values of A: A[,1] in the range <
  # -800 and below while other values in A[,1] were at least -400 and
  # mostly in range -100 ... 0; for A[,2] typical range around -800,
  # and for the NaN outlier: -1000. In both cases it was A[209,1] and
  # A[209,2] i.e. the same sample.  This is so rare that the effect on
  # the results is negligible, but this should be diagnosed and fixed
  # asap. As this is part of iteration, simply replace the defected sample with
  # equal probability in all groups.
  if (sum(!is.na(qOFz)) > 0) {
    # FIXME: optimize with apply!
    # Detect empty components and ignore
    inds <- c()
    for (i in 1:ncol(qOFz)) {
      if (sum(na.omit(qOFz[, i])) == 0) {
        inds <- c(inds, i)
        qOFz[, i] <- 0
      }
    }

    for (i in 1:nrow(qOFz)) {
      if (sum(is.na(qOFz[i, ])) > 0) {
        inds2 <- setdiff(seq(ncol(qOFz)), inds)
        qOFz[i, inds2] <- rep.int(1/length(inds2), length(inds2))
      }
    }
  }
  
  qOFz
  
}

mk.hp.posterior <- function(data, qOFz, hp.prior, opts){

  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  dat <- data$given.data$X1
  
  # Ensure that qOFz is a matrix
  qOFz <- matrix(qOFz, nrow(dat))

  # If qOFz exceeds max cluster size then remove one cluster 
  # (the second last one, which assumes clusters are sorted and the last is a new one)
  # FIXME 3/2012 quick hack - consider better implementations regarding this
  if (ncol(qOFz) > opts$c.max + 1) {
    inds <- setdiff(1:ncol(qOFz), ncol(qOFz) - 1)
    qOFz <- matrix(qOFz[, inds], nrow(dat))
    qOFz <- matrix(qOFz/rowSums(qOFz), nrow(dat))
  }

  # Compatibility variables not needed for the current functionality
  tmp.realS <- X2 <- dimX2 <- 0

  out <- .Call("mHPpost",
               dat,
               ncol(dat),
               nrow(dat),
               X2, dimX2,
               tmp.realS, opts$implicitnoisevar,
               hp.prior$Mumu, hp.prior$S2mu,
               hp.prior$AlphaKsi, hp.prior$BetaKsi,
               hp.prior$U.p, hp.prior$alpha,
               qOFz, ncol(qOFz), PACKAGE = "netresponse")

  qOFz <- matrix(out$qOFz, nrow(qOFz))

  #if (ncol(qOFz) > opts$c.max) {
  #}

  hp.posterior <- list(
    Mubar     = matrix(out$Mubar,    ncol(qOFz)),
    Mutilde   = matrix(out$Mutilde,  ncol(qOFz)),
    KsiAlpha  = matrix(out$KsiAlpha, ncol(qOFz)),
    KsiBeta   = matrix(out$KsiBeta,  ncol(qOFz)),
    gamma     = matrix(out$gamma, 2),
    Nc        = out$Nc,
    qOFz      = qOFz,
    Uhat      = out$Uhat)

  hp.posterior
}


mk.hp.prior <- function(data, opts){

  #  Copyright (C) 2008-2012 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

 # INPUT:   data     - matrix with data vectors
 #          opts     - list with algorithm options
 # OUTPUT:  hp_prior - list with prior information
 #
 # DESCRIPTION: Prior for the mean     = mean of data
 #              Prior for the variance = variance of data


  dat <- data$given.data$X1 # real-valued. Data to be clustered.

  Mean  <- colMeans(dat)     # mean of each dimension
  Var   <- colVariances(dat, Mean) # Variance of each dimension
                  #colSums((dat - rep(Mean, each = nrow(dat)))^2)/nrow(dat) 

  # priors for distribution of codebook vectors Mu ~ N(MuMu, S2.Mu)..
  #list(Mumu = Mean, S2mu = Var, U.p = Inf)
  # priors for data variance Ksi ~ invgam(AlphaKsi, BetaKsi)
  # variance is modeled with inverse Gamma distribution
  # FIXME: some of these are redundant, remove to save memory
  list(Mumu = Mean, S2mu = Var, U.p = Inf, AlphaKsi = rep(opts$prior.alphaKsi, ncol(dat)), BetaKsi = rep(opts$prior.betaKsi, ncol(dat)), alpha = opts$prior.alpha)
  
}


########################################################################################


# INPUT:   data, hp_posterior, hp_prior, opts
# OUTPUT:  free_energy: value of mixture model's free energy
#          log_lambda: Used for posterior of labels q_of_z <- softmax(log_lambda);
# DESCRIPTION: ...


mk.free.energy <- function(data, hp.posterior, hp.prior, opts,
                           fc = NULL, log.lambda = NULL)
{

  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner
  
  if( is.null(fc) || is.null(log.lambda) ){
    fc <- mk.E.log.q.p.eta(data, hp.posterior, hp.prior, opts)      # 1*K
    log.lambda <- mk.log.lambda(data, hp.posterior, hp.prior, opts) # N*K
  } 

  hpgsum <- colSums(hp.posterior$gamma)
  dig <- digamma(hpgsum)
  
  E.log.p.of.V <- lgamma(hpgsum) -
      lgamma(1 + hp.prior$alpha) -
      colSums(lgamma(hp.posterior$gamma)) +
      lgamma(hp.prior$alpha) +
      ( (hp.posterior$gamma[1, ] - 1) * 
        (digamma(hp.posterior$gamma[1, ]) - dig)) +
      ( (hp.posterior$gamma[2, ] - hp.prior$alpha) *
        (digamma(hp.posterior$gamma[2, ]) - dig))

  extra.term <- sum(E.log.p.of.V)
  free.energy <- extra.term + sum(fc) - sumlogsumexp(log.lambda)

  # Return
  list(free.energy = free.energy, log.lambda = log.lambda)

}
  
#mk.free.energy.c <- cmpfun(mk.free.energy)

mk.qOFz <- function(data, hp.posterior, hp.prior, opts, log.lambda = NULL){

  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  if( is.null(log.lambda) ){
    log.lambda <- mk.log.lambda(data, hp.posterior, hp.prior, opts)
  }

  qOFz <- softmax( log.lambda )

  # Do not allow empty clusters
  as.matrix(qOFz[, !colSums(qOFz) == 0], nrow(log.lambda))
  
}



mk.log.lambda <- function(data, hp.posterior, hp.prior, opts){


  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner

  dat <- data$given.data$X1
  
  # Ensure matrix for data
  if (!is.matrix(dat)) {stop("Error in mk.log.lambda: dat should be a matrix!")}
  
  # Compatibility variables, not needed for current functionality
  tmp.realS <- X2 <- dimX2 <- 0

  out <- .Call("mLogLambda",
               dat,
               ncol(dat),
               nrow(dat),
               X2, dimX2,
               tmp.realS, opts$implicitnoisevar,
               hp.prior, hp.posterior,
               PACKAGE = "netresponse")
  
  matrix(out, nrow(dat))
  
}


###########################################################################


# INPUT:   data, hp_posterior, hp_prior, opts
# OUTPUT:  matrix [1xk]: used to compute the free_energy formula. 
# DESCRIPTION: Regards the gaussian model's parameters.



mk.E.log.q.p.eta <- function(data, hp.posterior, hp.prior, opts){

  #  Copyright (C) 2008-2010 Antonio Gusmao and Leo Lahti
  #  Licence: GPL >=2
  #  This function is based on the Variational Dirichlet Process Gaussian
  #  Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
  #  (all rights reserved) and the Agglomerative Independent Variable
  #  Group Analysis package: Copyright (C) 2001-2007 Esa Alhoniemi,
  #  Antti Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and
  #  Paul Wagner
  
  # returns E [ log q(eta)/p(eta) ].q
  # fc: 1 by k

  dat <- data$given.data$X1
  
  # Ensure matrix for data
  if (!is.matrix(dat)) {stop("Error in mk.E.log.q.p.eta: dat is not a matrix!")}
  
  N  <- nrow(dat)
  M1 <- ncol(dat)
  K  <- nrow(hp.posterior$Mubar)
  l.codebook <- (- M1/2) * matrix(1, K)
  Ksi.log <- (digamma(hp.posterior$KsiAlpha) - log(hp.posterior$KsiBeta))
  
  for (j in 1:M1) {
    l.codebook <- l.codebook +
      .5 * ( log(hp.prior$S2mu[[j]] / hp.posterior$Mutilde[, j]) +
             ( (hp.posterior$Mubar[, j] - hp.prior$Mumu[[j]])^2 +
                hp.posterior$Mutilde[, j] 
             ) / hp.prior$S2mu[[j]]
           ) +
      lgamma(hp.prior$AlphaKsi[[j]]) - 
      lgamma(hp.posterior$KsiAlpha[, j]) +
      hp.posterior$KsiAlpha[, j] * log(hp.posterior$KsiBeta[, j]) -
      hp.prior$AlphaKsi[[j]] * log(hp.prior$BetaKsi[[j]]) + 
      (hp.posterior$KsiAlpha[, j] - hp.prior$AlphaKsi[[j]]) * Ksi.log[, j] +
      (hp.prior$BetaKsi[[j]] - hp.posterior$KsiBeta[, j]) *
        (hp.posterior$KsiAlpha[, j] / hp.posterior$KsiBeta[, j])
  }

  t(l.codebook)
}

#mk.E.log.q.p.eta.c <- cmpfun(mk.E.log.q.p.eta)

#################################################


colVariances <- function (dat, Mean) {
  # This is about 5x faster than apply(dat, 2, var)
  #max(abs(colVariances(dat) - apply(dat,2,var)))
  colSums((dat - rep(Mean, each = nrow(dat)))^2)/(nrow(dat) - 1)
}

#colVariances.c <- cmpfun(colVariances)

############################################################################

# INPUT:   data, qOFz, hp_posterior, hp_prior, opts
# OUTPUT:  list(new.qOFz, new.c);
#             * new.qOFz: posterior over labels including the split clusters.
#             * new.c: index of the newly created cluster.
# DESCRIPTION: Implements the VDP algorithm step 3a.

split.qofz <- function(qOFz, c, new.c, dat, speedup = TRUE, min.size = 4){

  # compute the first principal component of the candidate cluster,
  # not the whole data set.

  # min.size option sets the required minimum size of a cluster for
  # splitting; smaller clusters are not splitted

  # Pick sample indices and samples corresponding to cluster c
  cluster_assignments <- apply(qOFz, 1, which.max);
  indices <- which(cluster_assignments == c);
  if (length(indices) < min.size) {
    #"Component must have at least min.size samples to be splitted."
    # -> no splitting
    new.qOFz <- qOFz
  } else {

    component.data <- matrix(dat[indices,], length(indices))
  
    # If the number of samples is high calculating PCA might take long
    # but can be approximated by using less samples:

    pcadata <- component.data
    
    if ( speedup ) {

      # when a candidate cluster, C, is split to generate two new
      # clusters, it is split by mapping the data onto the first
      # principal component of the data in C and then splitting that in
      # half. To speed up, one can compute an approximate first
      # principal component by considering a randomly selected subset of
      # the data belonging to C, and computing its first principal
      # component.

      # number of samples in this component
      ns <- nrow(component.data) 
      nd <- ncol(component.data) 

      # If component size exceeds cmax, 
      # use only a random subset of data to calculate PCA
      # size of the random subset increases slowly (linearly) 
      # with component size. 
      cmax <- 20 #take at least this many random samples    
      nr <- min(ns, cmax + floor(sqrt(ns))) # do not take more samples than are available
      rinds <- sample(ns, nr)

      # Pick random subset of the component data and accompanying indices
      # to speed up PCA calculations
      pcadata <- matrix(component.data[rinds,], nrow = nr)
      indices <- indices[rinds]
    }

    # Split the cluster based on the first PCA component
    # FIXME: compare speed with other PCA implementations and select fastest
    dir <- prcomp(pcadata)$x[,1]
    I1 <- indices[dir >= 0];
    I2 <- indices[dir < 0];

    # Initialize split by adding a zero column on qOFz

    # If one of qOFz clusters is empty, then do not create new clusters but instead fill in the empty cluster
    # during cluster split.
    # FIXME: ensure already in creating qOFz-matrices that no zero columns are allowed. This will
    # avoid the need to address the issue here.
    # -> OK, done this. w remove this unnecessary check here and test if
    # the code works ok
    empty.cols <- (colSums(qOFz) == 0)
    if ( !any(empty.cols) ) { # no empty columns -> add an empty cluster
      new.qOFz <- array(0, dim = c(nrow(qOFz), ncol(qOFz) + 1))
      new.qOFz[,  -new.c] <- qOFz
    } else { # an empty column -> no need to add new clusters
      new.qOFz <- qOFz
      new.c <- which(empty.cols)[[1]]
    }

    # Split this component (samples given in I1, I2) into two smaller components
    new.qOFz[ I1, c]     <- qOFz[ I1, c]
    new.qOFz[ I2, c]     <- 0 # Remove entries from cluster c
    new.qOFz[ I2, new.c] <- qOFz[ I2, c] # Add same entries to cluster new.c
  }

  new.qOFz

}

#split.qofz.c <- cmpfun(split.qofz)

###############################################################################


#fsort <- function (df, sortvar, decreasing = FALSE) {
fsort <- function (df, sortvar) {

  o <- order(df[[sortvar]])
  #if (decreasing) {o <- rev(o)}
  df[o,]

}

#fsort.c <- cmpfun(fsort)

##################################################################
	

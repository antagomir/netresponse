get.model.parameters <- function (model, subnet.id, datamatrix, level = NULL) {
            
  # model: output from run.netresponse function
  # subnet.idx: index of the subnet to check
  # level: which agglomeration step
  
  #  Copyright (C) 2008-2010 Leo Lahti
  #  Licence: GPL >=2
 
  nodes <- get.subnets(model, level)[[subnet.id]]

  # Compute the model
  m <- vdp.mixt(matrix(datamatrix[, nodes], nrow = length(model@samples) ))

  # Pick parameters
  sds  <- sqrt(m$variances)              # component standard devs
  mu   <- m$means #m$hp.posterior$Mubar  # component centroids
  w    <- m$weights                      # component weights

  rownames(mu) <- rownames(sds) <- names(w) <- paste("Response", 1:length(w), sep = "-")
  colnames(mu) <- colnames(sds) <- nodes

  # For mu and std, rows correspond to the clusters, in w the elements
  list(mu = mu, std = sds, w = w, nodes = nodes, K = m$K)

}

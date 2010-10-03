get.model.parameters <- function (model, subnet.id, level = NULL) {
            
  # model: output from run.netresponse function
  # subnet.id: index of the subnet to check
  # level: which agglomeration step
  
  #  Copyright (C) 2008-2010 Leo Lahti
  #  Licence: GPL >=2
 
  nodes <- get.subnets(model, level)[[subnet.id]]
  datamatrix <- model@datamatrix  

  # Compute the model
  m <- vdp.mixt(matrix(datamatrix[, nodes], nrow = length(model@samples) ))

  # Pick parameters
  w    <- m$posterior$weights    # component weights
  mu   <- m$posterior$centroids  # component centroids
  sds  <- sqrt(m$posterior$sds)  # component standard devs

  rownames(mu) <- rownames(sds) <- names(w) <- paste("Response", 1:length(w), sep = "-")
  colnames(mu) <- colnames(sds) <- nodes

  # For mu and std, rows correspond to the mixture components, in w the elements
  list(mu = mu, sd = sds, w = w, nodes = nodes, K = m$K)

}

model.stats <- function ( model ) {

  # Check statistics for subnetworks
  # subnetwork size
  # number of responses

  subnets <- model@subnets

  Ncomps <- c()
  for (subnet.id in names(subnets)) {
    #print(subnet.id)
    # number of mixture components
    Ncomps[[subnet.id]] <- model@models[[subnet.id]]$K
    
  }

  tab <- cbind(sapply(subnets, length), Ncomps)  
  colnames(tab) <- c("subnet.size", "subnet.responses")
  rownames(tab) <- paste("Subnet", 1:nrow(tab), sep="-")

  as.data.frame(tab)

}


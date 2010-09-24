result.stats <-
  
function ( model, datamatrix, level = NULL ) {


  # Check statistics for subnetworks
  # subnetwork size
  # number of responses

  # FIXME add level option
  subnets <- model@last.grouping
  
  Ncomps <- c()
  for (subnet.id in names(subnets)) {
    
    m <- get.model(model, subnet.id, datamatrix, level)
    
    # number of mixture components
    Ncomps[[subnet.id]] <- m$K
    
  }

  tab <- cbind(sapply(subnets, length), Ncomps)  
  colnames(tab) <- c("subnet.size", "subnet.responses")
  rownames(tab) <- paste("Subnet", 1:nrow(tab), sep="-")

  as.data.frame(tab)

}


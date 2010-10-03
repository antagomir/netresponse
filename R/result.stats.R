result.stats <-
  
function ( model, level = NULL ) {


  # Check statistics for subnetworks
  # subnetwork size
  # number of responses

  # FIXME add level option
  subnets <- model@last.grouping
  datamatrix <- model@datamatrix  

  Ncomps <- c()
  for (subnet.id in names(subnets)) {
    
    m <- get.model(model, subnet.id, level)
    
    # number of mixture components
    Ncomps[[subnet.id]] <- m$posterior$K
    
  }

  tab <- cbind(sapply(subnets, length), Ncomps)  
  colnames(tab) <- c("subnet.size", "subnet.responses")
  rownames(tab) <- paste("Subnet", 1:nrow(tab), sep="-")

  as.data.frame(tab)

}


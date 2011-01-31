response2sample <-
function (model, subnet.id, component.list = TRUE) {

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", "Subnet-", subnet.id, sep="")
  }
  
  response.probabilities <- sample2response(model, subnet.id)

  # For each sample, list the most strongly associated response (highest P(r|s))
  clusters <- apply(response.probabilities, 1, which.max)
  
  if ( component.list ) {
    # list samples separately for each cluster
    clusters <- lapply(seq(max(clusters)), function( i ){ names(which(clusters == i)) })
    # names(clusters) <- FIXME add names here
    if (length(clusters) < ncol(response.probabilities)) {
      n <- ncol(response.probabilities) - length(clusters)
      clusters <- c(clusters, vector(n, mode = "list"))
    }
  }
    
  clusters
  
}


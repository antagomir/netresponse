response2sample <-
function (model, subnet.id, datamatrix, component.list = TRUE) {
  response.probabilities <- sample2response(model, subnet.id, datamatrix)

  # For each sample, list the most strongly associated response (highest P(r|s))
  clusters <- apply(response.probabilities, 1, which.max)
  
  if ( component.list ) {
    # list samples separately for each cluster
    clusters <- lapply(seq(max(clusters)), function( i ){ names(which(clusters == i)) })
  }
  
  clusters
  
}


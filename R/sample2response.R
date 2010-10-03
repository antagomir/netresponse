sample2response <-
function (model, subnet.id) {

  # Retrieve model for given subnet
  m <- get.model(model, subnet.id) 

  # P(response | sample)
  assignment.matrix <- m$posterior$qOFz
  rownames(assignment.matrix) <- model@samples
  
  assignment.matrix
  
}


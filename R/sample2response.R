sample2response <-
function (model, subnet.id, datamatrix) {

  # Retrieve model for given subnet
  m <- get.model(model, subnet.id, datamatrix) 

  # P(response | sample)
  assignment.matrix <- m$qOFz
  rownames(assignment.matrix) <- model@samples
  
  assignment.matrix
  
}


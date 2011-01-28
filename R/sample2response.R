sample2response <- function (model, subnet.id) {

  # P(response | sample)
  assignment.matrix <- model@models[[subnet.id]]$posterior$qOFz
  rownames(assignment.matrix) <- model@samples
  colnames(assignment.matrix) <- paste("Response", 1:ncol(assignment.matrix), sep = "-")
  
  assignment.matrix
  
}


sample2response <- function (model, subnet.id) {

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", "Subnet-", subnet.id, sep="")
  }
  
  # P(response | sample)
  assignment.matrix <- model@models[[subnet.id]]$posterior$qOFz
  rownames(assignment.matrix) <- model@samples
  colnames(assignment.matrix) <- paste("Response", 1:ncol(assignment.matrix), sep = "-")
  
  assignment.matrix
  
}


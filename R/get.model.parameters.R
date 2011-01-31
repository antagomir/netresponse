get.model.parameters <- function (model, subnet.id) {
            
  # model: output from run.netresponse function
  # subnet.id: index of the subnet to check
  
  #  Copyright (C) 2008-2011 Leo Lahti
  #  Licence: GPL >=2

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", "Subnet-", subnet.id, sep="")
  }
  
  pars <- model@models[[subnet.id]]
  pars[["nodes"]] <- model@subnets[[subnet.id]]

  pars

}

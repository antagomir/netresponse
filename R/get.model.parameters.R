get.model.parameters <- function (model, subnet.id) {
            
  # model: output from run.netresponse function
  # subnet.id: index of the subnet to check
  
  #  Copyright (C) 2008-2011 Leo Lahti
  #  Licence: GPL >=2
 
  pars <- model@models[[subnet.id]]
  pars[["nodes"]] <- model@subnets[[subnet.id]]

  pars

}

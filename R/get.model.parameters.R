
# Copyright (C) 2010-2012 Leo Lahti
# Contact: Leo Lahti <leo.lahti@iki.fi>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

get.model.parameters <- function (model, subnet.id) {
            
  # model: output from run.netresponse function
  # subnet.id: index of the subnet to check
  
  #  Copyright (C) 2008-2011 Leo Lahti
  #  Licence: GPL >=2

  if (is.numeric(subnet.id)) {
    warning("subnet.id given as numeric; converting to character: ", "Subnet-", subnet.id, sep="")    
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
  }
  
  pars <- model@models[[subnet.id]]
  pars[["nodes"]] <- model@subnets[[subnet.id]]

  pars

}

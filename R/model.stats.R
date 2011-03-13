
# Copyright (C) 2010-2011 Leo Lahti
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



model.stats <- function ( model ) {

  # Check statistics for subnetworks
  # subnetwork size
  # number of responses

  subnets <- model@subnets

  Ncomps <- c()
  for (subnet.id in names(subnets)) {
    #print(subnet.id)
    # number of mixture components
    Ncomps[[subnet.id]] <- length(model@models[[subnet.id]]$w)
  }

  tab <- cbind(sapply(subnets, length), Ncomps)  
  colnames(tab) <- c("subnet.size", "subnet.responses")
  rownames(tab) <- names(subnets)

  as.data.frame(tab)

}


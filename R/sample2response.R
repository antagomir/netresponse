
# Copyright (C) 2010-2012 Leo Lahti
# Contact: Leo Lahti <leo.lahti@iki.fi>
# This file is part of NetResponse R program
#  
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

sample2response <- function (model, subnet.id, mode = "soft") {

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", subnet.id, sep="")
  }
  
  # P(response | sample)
  #assignment.matrix <- model@models[[subnet.id]]$posterior$qOFz
  assignment.matrix <- getqofz(model, subnet.id, log = FALSE)
  
  if (mode == "hard") {
    sample.names <- rownames(assignment.matrix)
    assignment.matrix <- colnames(assignment.matrix)[apply(assignment.matrix, 1, which.max)]
    names(assignment.matrix) <- sample.names
  }

  assignment.matrix

}



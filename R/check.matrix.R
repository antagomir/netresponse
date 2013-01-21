# Copyright (C) 2008-2012 Olli-Pekka Huovilainen and Leo Lahti 
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

# The best academic advice I ever got was: "Spend at least an hour
# every day on the manuscript closest to publication".

#' check.matrix
#' 
#' Mostly for internal purposes. Check input matrix format.
#' 
#' 
#' @usage check.matrix(datamatrix)
#' @param datamatrix See detect.responses
#' @return The datamatrix, possibly added with necessary formatting for the
#' netresponse algorithm.
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.  Maintainer:
#' Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso detect.responses
#' @references See citation("netresponse")
#' @keywords internal
#' @examples # datamatrix <- check.matrix(datamatrix)
check.matrix <- function (datamatrix) {
  
  accepted.formats.emat <- c("matrix", "Matrix", "data.frame")  

  # ensure datamatrix is a matrix
  if (!is.matrix(datamatrix)) {
    if (class(datamatrix) %in% accepted.formats.emat) {
      datamatrix <- as.matrix(datamatrix)
    } else {
      stop(paste("datamatrix needs to be in one of the following formats:", paste(accepted.formats.emat, collapse = "; ")))
    }    
  }
  if (is.null(colnames(datamatrix))) { colnames(datamatrix) <- as.character(1:ncol(datamatrix)) }
  if (is.null(rownames(datamatrix))) { rownames(datamatrix) <- as.character(1:nrow(datamatrix)) }  

  datamatrix
}


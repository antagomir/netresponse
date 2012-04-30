#' check.matrix
#' 
#' Mostl for internal purposes. Check input matrix format.
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


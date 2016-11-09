#' @title BIC mixture
#' @description Latent class analysis based on (infinite) Gaussian mixture model. If the input is data matrix, a multivariate model is fitted; if the input is a vector, a univariate model is fitted
#' @param x samples x features matrix for multivariate analysis, or a vector for univariate analysis 
#' @param max.modes Maximum number of modes to be checked for mixture model selection
#' @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture.
#' @param min.modes minimum number of modes
#' @param ... Further optional arguments to be passed
#' @return Fitted latent class model (parameters and free energy)
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
bic.mixture <- function (x, max.modes, bic.threshold = 0, min.modes = 1, ...) { 

   # x; max.modes = max.responses; bic.threshold = bic.threshold; min.modes = min.responses

  if (!is.vector(x) && ncol(x) == 1) {x <- x[,1]}	    

  if (is.vector(x)) {
    bic.mixture.univariate(x, max.modes, bic.threshold, min.modes = min.modes, ...)
  } else {
    bic.mixture.multivariate(x, max.modes, bic.threshold, min.modes = min.modes, ...)
  }

}



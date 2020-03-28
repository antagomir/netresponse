#' @title Add ellipse to an existing plot
#' @description Calculates and plots ellipse corresponding to specified confidence interval in 2-dimensional plot
#' @usage add.ellipse(centroid, covmat, confidence = 0.95, npoints = 100, col =
#' "black", ...)
#' @param centroid Vector with two elements defining the ellipse centroid.
#' @param covmat Covariance matrix for the investigated data. Only diagonal
#' covariances supported.
#' @param confidence Confidence level determining the ellipse borders based on
#' the covariance matrix.
#' @param npoints Number of plotting points.
#' @param col Color.
#' @param ... Other arguments to be passed.
#' @return Used for plotting side effects.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
#' @export
#' @examples #add.ellipse(centroid = c(0, 0), covmat = diag(c(1,2)))
add.ellipse <- function (centroid, covmat, confidence = 0.95, npoints = 100, col = "black", ...) {

  # add ellipse to a plot 
  el <- ellipse(centroid, covmat, confidence, npoints)
  points(el, type = "l", col = col, ...)

  el
}




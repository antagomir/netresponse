# "The language of science is the language of probability, and not of
#  p-values." -- Luis Pericchi

#' @title Set breaks
#' @description Set breakpoints for two-way color palette.
#' @usage set.breaks(mat, interval = 0.1)
#' @param mat Matrix to visualize.
#' @param interval Density of color breakpoints.
#' @return A vector listing the color breakpoints.
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao. Maintainer:
#' Leo Lahti <leo.lahti@@iki.fi>
#' @references L. Lahti et al.: Global modeling of transcriptional responses in
#' interaction networks. Submitted.
#' @keywords utilities
#' @export
#' @examples
#'  set.breaks(array(rnorm(100), dim = c(10, 10)), interval = .1) 
set.breaks <- function (mat, interval = .1) {
  if (max(abs(mat)) > 1) {
    m <- floor(max(abs(mat)))
  } else {
    m <- round(max(abs(mat)), nchar(1/interval) - 1)
  }

  mm <- m + interval/2
  vals <- seq(interval/2,mm,interval)
  # Note: the first and last values mimic infinity
  mybreaks  <- c(-(m + 1e6), c(-rev( vals ), vals), m + 1e6)
  mybreaks
}



ellipse <- function (centroid, covmat, confidence = 0.95, npoints = 100) {

  # Centroid: x0, y0
  # axes: a, b
  x0 <- centroid[[1]]
  y0 <- centroid[[2]]

  # Pick axis-wise stds
  a <- sqrt(diag(covmat))[[1]]
  b <- sqrt(diag(covmat))[[2]]

  theta <- seq(0, 2 * pi, length = npoints)
 
  # Confidence intervals (df=1: 1.96; df=2: 2.45)
  cint <- sqrt(qchisq(confidence, 2))

  # Determine point on the ellipse
  x <- x0 + cint * a * cos(theta)
  y <- y0 + cint * b * sin(theta)

  # Rotate alpha degrees from the x-axis:
  #x <- x0 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
  #y <- y0 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)

  # Output ellipse points
  cbind(x, y)

}




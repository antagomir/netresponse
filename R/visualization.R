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


#' PlotMixtureUnivariate
#' 
#' Visualize data, centroids and stds for a given univariate
#' Gaussian mixture model with PCA.
#'
#' Arguments:
#' @param x data vector
#' @param means mode centroids
#' @param sds mode standard deviations
#' @param ws weight for each mode
#' @param title.text Plot title
#' @param xlab.text xlab.text
#' @param ylab.text ylab.text
#' @param binwidth binwidth for histogram
#' @param ... Further arguments for plot function.
#'
#' Return:
#' @return Used for its side-effects
#'
#' @export
#'
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse") for citation details.
#' @keywords utilities
#' @examples #plotMixtureUnivariate(dat, means, sds, ws)
PlotMixtureUnivariate <- function (x, means, sds, ws, title.text = NULL, xlab.text = NULL, ylab.text = NULL, binwidth = 0.05, ...) {

  # Circumvent warnings
  ..density.. <- NULL
  vals <- NULL
  varname <- NULL
		      
  df <- data.frame(list(x = x))
  df$vals <- seq(min(df$x), max(df$x), length=nrow(df)) # estimation points for fitted Gaussians

  # Histogram and density plot
  pg <- ggplot2::ggplot(df, aes(x=x)) 
  pg <- pg + geom_histogram(aes(y = ..density..), binwidth=binwidth, fill = "gray") 
  pg <- pg + geom_density(fill="gray", alpha = 0.1) 
  pg <- pg + theme_bw() + xlab(xlab.text) + ylab(ylab.text) 
  pg <- pg + ggtitle(title.text)

  # Estimated normal distributions from the mixture model
  for(comp in 1:length(means)){
    df2 <- data.frame(list(vals = df$vals, varname = ws[[comp]]*dnorm(df$vals, mean = means[[comp]], sd = sds[[comp]])))
    pg <- pg + geom_line(aes(x=vals, y=varname), colour="red", data = df2) # Add in normal distribution
  }

  pg

}




# "The language of science is the language of probability, and not of
#  p-values." -- Luis Pericchi

#' set.breaks
#' 
#' Set breakpoints for two-way color palette.
#' 
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



#' Add ellipse to an existing plot.
#' 
#' Calculates and plots ellipse corresponding to specified confidence interval
#' in 2-dimensional plot.
#' 
#' 
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
#' @examples
#' 
#' #add.ellipse(centroid = c(0, 0), covmat = diag(c(1,2)))
#' 
add.ellipse <- function (centroid, covmat, confidence = 0.95, npoints = 100, col = "black", ...) {

  # add ellipse to a plot 
  el <- ellipse(centroid, covmat, confidence, npoints)
  points(el, type = "l", col = col, ...)

  el
}



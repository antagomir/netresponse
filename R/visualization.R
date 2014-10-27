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

plotMatrix.2way <- function (mat, mybreaks = NULL, maintext = "", xlab = "", ylab = "", mypalette = NULL, interval = .1, cex.main = 1, xaxis = FALSE, yaxis = TRUE, row.tick = 1, col.tick = 1, cex.xlab = .9, cex.ylab = .9, cex.lab = .9, limit.trunc = 0, mar = c(5, 4, 4, 2), ...) {

  # mat: differential expression matrix to plot in two-color palette
  # interval: interval for palette color switches
  # FIXME: synchronize with PlotMatrix in sorvi package  
   	   
  if (length(mybreaks) == 0)  {
    m <- max(round(max(abs(mat)), limit.trunc) - interval, 0)
    mm <- m + interval/2
    vals <- seq(interval/2,mm,interval)
    # Set breaks evenly around zero
    mybreaks  <- c(-(m+1e6),c(-rev(vals),vals),m+1e6)
  }
		  
  if (length(mypalette)==0) {
    mypalette <- colorRampPalette(c("blue", "black", "red"),space = "rgb")
    my.colors <- mypalette(length(mybreaks)-1)
  } else {
    my.colors <- mypalette(length(mybreaks)-1)
  }
		      
  # transpose and revert row order to plot matrix in the same way it
  # appears in its numeric form
  par(mar = mar)
  image(t(mat[rev(seq(nrow(mat))),]), col = my.colors, xaxt='n', yaxt='n', zlim=range(mybreaks), breaks=mybreaks, main=maintext, xlab=xlab, ylab=ylab, cex.lab = cex.lab, cex.main = cex.main)

  if (yaxis) {
      
    v <- seq(1, nrow(mat), row.tick) # take every nth index
    axis(2, at = seq(0,1,length = nrow(mat))[v], labels = rev(rownames(mat))[v], 
    	    cex.axis=cex.ylab, las = 2)
    
  }
  
  if (xaxis) {    

    v <- seq(1, ncol(mat), col.tick) # take every nth index
    axis(1, at = seq(0,1,length = ncol(mat))[v], labels = colnames(mat)[v], 
    	    cex.axis = cex.xlab, las=2)

  }
    
  return(list(palette = my.colors, breaks = mybreaks))
      	  
}


check.bins <- function (difexp, mybreaks) {

  # check color scale bin for each expression value
  bins <- c()
  for (i in 1:length(difexp)) {
    # which color bins are smaller than our difexp value
    # (for probet: i, mode:mode)
    inds <- which(difexp[[i]] > mybreaks)
    if (length(inds) == 0) {
      bins[[i]] <- 1
    } else if (length(inds) > 0)  {
      bins[[i]] <- max(inds) + 1
    }
  }	
  bins
}


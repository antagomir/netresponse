
write.netresponse.results <- function (x, subnet.ids = NULL, filename) {

  #f <- file(description = filename, open = "rw")   
  write("NetResponse - subnetworks", file = filename, append = FALSE)
  write("==========================\n", file = filename, append = TRUE)

  if ( is.null(subnet.ids) ) { subnet.ids <- names(x@subnets) }

  for (nam in subnet.ids) {
    write(nam, file = filename, append = TRUE)    
    write(x@subnets[[nam]], file = filename, append = TRUE)
    write("----------------", file = filename, append = TRUE)    
  }

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

add.ellipse <- function (centroid, covmat, confidence = 0.95, npoints = 100, col = "black", ...) {

  # add ellipse to a plot 
  el <- ellipse(centroid, covmat, confidence, npoints)
  points(el, type = "l", col = col, ...)

  el
}


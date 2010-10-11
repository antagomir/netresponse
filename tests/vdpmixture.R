
# 1. vdp.mixt: moodien loytyminen eri dimensiolla, naytemaarilla ja komponenteilla
#   -> ainakin nopea check

#######################################################################

# Generate random data from five Gaussians. 
# Detect modes with vdp-gm. 
# Plot data points and detected clusters with variance ellipses

#######################################################################

library(netresponse)

#dyn.load("/home/tuli/Rpackages/netresponse/netresponse/src/netresponse.so")


###################################################################



ellipse <-
  function (x, scale = c(1, 1), centre = c(0, 0), level = 0.95, 
            t = sqrt(qchisq(level, 2)), which = c(1, 2), npoints = 100, ...) 
{

  # Modified from the 'ellipse.default' function of ellipse package (version 0.3-5) in CRAN

  names <- c("x", "y")
  if (is.matrix(x)) {
    xind <- which[1]
    yind <- which[2]
    r <- x[xind, yind]
    if (missing(scale)) {
      scale <- sqrt(c(x[xind, xind], x[yind, yind]))
      if (scale[1] > 0) r <- r/scale[1]
      if (scale[2] > 0) r <- r/scale[2]
    }
    if (!is.null(dimnames(x)[[1]])) 
      names <- dimnames(x)[[1]][c(xind, yind)]
  }
  else r <- x
  r <- min(max(r,-1),1)  # clamp to -1..1, in case of rounding errors
  d <- acos(r)
  a <- seq(0, 2 * pi, len = npoints)
  matrix(c(t * scale[1] * cos(a + d/2) + centre[1], t * scale[2] * 
           cos(a - d/2) + centre[2]), npoints, 2, dimnames = list(NULL, 
                                                    names))
}



#########  GENERATE DATA #############################################

# Generate Nc components from normal-inverseGamma prior

set.seed(12346)

dd <- 2   # Dimensionality of data
Nc <- 5   # Number of components
Ns <- 200 # Number of data points
sd0 <- 3  # component spread
rgam.shape = 2 # parameters for Gamma distribution 
rgam.scale = 2 # parameters for Gamma distribution to define precisions


# Generate means and variances (covariance diagonals) for the components 
component.means <- matrix(rnorm(Nc*dd, mean = 0, sd = sd0), nrow = Nc, ncol = dd)
component.vars <- matrix(1/rgamma(Nc*dd, shape = rgam.shape, scale = rgam.scale), 
	                 nrow = Nc, ncol = dd)
component.sds <- sqrt(component.vars)


# Size for each component -> sample randomly for each data point from uniform distr.
# i.e. cluster assignments
sample2comp <- sample.int(Nc, Ns, replace = TRUE)

D <- array(NA, dim = c(Ns, dd))
for (i in 1:Ns)  {
    # component identity of this sample
    ci <- sample2comp[[i]]
    cm <- component.means[ci,]
    csd <- component.sds[ci,]
    D[i,] <- rnorm(dd, mean = cm, sd = csd)
}


######################################################################

# Fit nonparametric Gaussian mixture model
out <- vdp.mixt(D)

############################################################

# Compare input data and results

ord.out <- order(out$posterior$centroids[,1])
ord.in <- order(component.means[,1])

means.out <- out$posterior$centroids[ord.out,]
means.in <- component.means[ord.in,]

# Cluster stds and variances
sds.out <- out$posterior$sds[ord.out,]
sds.in  <- sqrt(component.vars[ord.in,])
vars.out <- sds.out^2
vars.in <- sds.in^2

# Check correspondence between input and output
if (length(means.in) == length(means.out)) {
   cm <- cor(as.vector(means.in), as.vector(means.out))
   csd <- cor(as.vector(sds.in), as.vector(sds.out))
}

# Plot results (assuming 2D)

ran <- range(c(as.vector(means.in - 2*vars.in), 
               as.vector(means.in + 2*vars.in), 
	       as.vector(means.out + 2*vars.out), 
	       as.vector(means.out - 2*vars.out)))

plot(D, pch = 20, main = paste("Cor.means:", round(cm,3), "/ Cor.sds:", round(csd,3)), xlim = ran, ylim = ran) 

for (ci in 1:nrow(means.out))  {
    points(means.out[ci,1], means.out[ci,2], col = "red", pch = 19)
    el <- ellipse(matrix(c(vars.out[ci,1],0,0,vars.out[ci,2]),2), centre = means.out[ci,])
    lines(el, col = "red") 						  
}

for (ci in 1:nrow(means.in))  {
    points(means.in[ci,1], means.in[ci,2], col = "blue", pch = 19)
    el <- ellipse(matrix(c(vars.in[ci,1],0,0,vars.in[ci,2]),2), centre = means.in[ci,])
    lines(el, col = "blue") 						  
}






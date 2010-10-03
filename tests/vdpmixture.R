
# 1. vdp.mixt: moodien loytyminen eri dimensiolla, naytemaarilla ja komponenteilla
#   -> ainakin nopea check

#######################################################################

# Generate random data from five Gaussians. 
# Detect modes with vdp-gm. 
# Plot data points and detected clusters with variance ellipses

#######################################################################

library(netresponse)

#dyn.load("/home/tuli/Rpackages/netresponse/netresponse/src/netresponse.so")

#########  GENERATE DATA #############################################

# Generate Nc components from normal-inverseGamma prior

set.seed(1234)

dd <- 2   # Dimensionality of data
Nc <- 3   # Number of components
Ns <- 100 # Number of data points
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

# Cluster variances
sds.out <- out$posterior$sds[ord.out,]
sds.in  <- sqrt(component.vars[ord.in,])

# Check correspondence between input and output
cm <- cor(as.vector(means.in), as.vector(means.out))
csd <- cor(as.vector(sds.in), as.vector(sds.out))

# Plot results (assuming 2D)
plot(D, pch = 20, main = paste("Cor.means:", round(cm,3), "/ Cor.sds:", round(csd,3))) 
for (ci in 1:Nc)  {
    points(means.out[ci,1], means.out[ci,2], col = "red", pch = 19)
    points(means.in[ci,1], means.in[ci,2], col = "blue", pch = 19)
}



# Generate Nc components from normal-inverseGamma prior

set.seed(12346)

dd <- 3   # Dimensionality of data
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


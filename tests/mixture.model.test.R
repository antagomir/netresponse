
# Validate mixture models


# Generate random data from five Gaussians. 
# Detect modes 
# Plot data points and detected clusters 

library(netresponse)

#fs <- list.files("~/Rpackages/netresponse/netresponse/R/", full.names = TRUE); for (f in fs) {source(f)}; dyn.load("/home/tuli/Rpackages/netresponse/netresponse/src/netresponse.so")

#########  Generate DATA #######################

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

######################################################################

par(mfrow = c(2,1))

for (mm in c("vdp", "bic")) {

# Fit nonparametric Gaussian mixture model
out <- mixture.model(D, mixture.method = mm, c.max = 10, pca.basis = FALSE)

############################################################

# Compare input data and results

ord.out <- order(out$params$mu[,1])
ord.in <- order(component.means[,1])

means.out <- out$params$mu[ord.out,]
means.in <- component.means[ord.in,]

# Cluster stds and variances
sds.out <- out$params$sd[ord.out,]
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

real.modes <- sample2comp
obs.modes <- apply(out$params$qofz, 1, which.max)

#plot(D, pch = 20, main = paste(mm, "/ cor.means:", round(cm,6), "/ Cor.sds:", round(csd,6)), xlim = ran, ylim = ran) 
plot(D, pch = real.modes, col = obs.modes, main = paste(mm, "/ cor.means:", round(cm,6), "/ Cor.sds:", round(csd,6)), xlim = ran, ylim = ran) 
for (ci in 1:nrow(means.out))  { add.ellipse(centroid = means.out[ci,], covmat = diag(vars.out[ci,]), col = "red") }
for (ci in 1:nrow(means.in))  { add.ellipse(centroid = means.in[ci,], covmat = diag(vars.in[ci,]), col = "blue") }

}

#######################################################

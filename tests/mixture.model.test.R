# Validate mixture models

# Generate random data from five Gaussians. 
# Detect modes 
# Plot data points and detected clusters 

library(netresponse)

#fs <- list.files("~/Rpackages/netresponse/netresponse/R/", full.names = TRUE); for (f in fs) {source(f)}; dyn.load("/home/tuli/Rpackages/netresponse/netresponse/src/netresponse.so")

#########  Generate DATA #######################

res <- generate.toydata()
D <- res$data
component.means <- res$means
component.sds   <- res$sds
sample2comp     <- res$sample2comp

######################################################################

par(mfrow = c(2,1))

for (mm in c("vdp", "bic")) {

  # Fit nonparametric Gaussian mixture model
  #source("~/Rpackages/netresponse/netresponse/R/vdp.mixt.R")
  out <- mixture.model(D, mixture.method = mm, max.responses = 10, pca.basis = FALSE)

  ############################################################

  # Compare input data and results

  ord.out <- order(out$mu[,1])
  ord.in <- order(component.means[,1])

  means.out <- out$mu[ord.out,]
  means.in <- component.means[ord.in,]

  # Cluster stds and variances
  sds.out <- out$sd[ord.out,]
  vars.out <- sds.out^2

  sds.in  <- component.sds[ord.in,]
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
  obs.modes <- apply(out$qofz, 1, which.max)

  # plot(D, pch = 20, main = paste(mm, "/ cor.means:", round(cm,6), "/ Cor.sds:", round(csd,6)), xlim = ran, ylim = ran) 
  plot(D, pch = real.modes, col = obs.modes, main = paste(mm, "/ cor.means:", round(cm,6), "/ Cor.sds:", round(csd,6)), xlim = ran, ylim = ran) 
  for (ci in 1:nrow(means.out))  { add.ellipse(centroid = means.out[ci,], covmat = diag(vars.out[ci,]), col = "red") }
  for (ci in 1:nrow(means.in))  { add.ellipse(centroid = means.in[ci,], covmat = diag(vars.in[ci,]), col = "blue") }

}


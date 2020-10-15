
# 1. vdp.mixt: moodien loytyminen eri dimensiolla, naytemaarilla ja komponenteilla
#   -> ainakin nopea check

#######################################################################

# Generate random data from five Gaussians. 
# Detect modes with vdp-gm. 
# Plot data points and detected clusters with variance ellipses

#######################################################################

library(netresponse)
#source("~/Rpackages/netresponse/netresponse/R/detect.responses.R")
#source("~/Rpackages/netresponse/netresponse/R/internals.R")
#source("~/Rpackages/netresponse/netresponse/R/vdp.mixt.R")
#dyn.load("/home/tuli/Rpackages/netresponse/netresponse/src/netresponse.so")


#########  Generate DATA #############################################

res <- generate.toydata()
D <- res$data
component.means <- res$means
component.sds   <- res$sds
sample2comp     <- res$sample2comp

######################################################################

# Fit nonparametric Gaussian mixture model
out <- vdp.mixt(D)
# out <- vdp.mixt(D, c.max = 3) # try with limited number of components -> OK

############################################################

# Compare input data and results

ord.out <- order(out$posterior$centroids[,1])
ord.in <- order(component.means[,1])

means.out <- out$posterior$centroids[ord.out,]
means.in <- component.means[ord.in,]

# Cluster stds and variances
sds.out <- out$posterior$sds[ord.out,]
sds.in  <- component.sds[ord.in,]
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

#plot(D, pch = 20, main = paste("Cor.means:", round(cm,3), "/ Cor.sds:", round(csd,3)), xlim = ran, ylim = ran) 
#for (ci in 1:nrow(means.out))  { add.ellipse(centroid = means.out[ci,], covmat = diag(vars.out[ci,]), col = "red") }
#for (ci in 1:nrow(means.in))  { add.ellipse(centroid = means.in[ci,], covmat = diag(vars.in[ci,]), col = "blue") }



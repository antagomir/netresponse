# Test script for the ICMg method

# Load the package
library(netresponse)

data(osmo) # Load data

# Set parameters
C.boost = 1
alpha = 10
beta = 0.01
B.num = 10
B.size = 10
S.num = 10  
S.size = 10
C = 24
pm0 = 0
V0 = 1               
V = 0.1

# Run combined ICMg sampler
res = ICMg.combined.sampler(osmo$ppi, osmo$exp, C, alpha, beta, pm0, V0, V, B.num, B.size, S.num, S.size, C.boost) 

# Compute component membership probabilities for the data points
res$comp.memb <- ICMg.get.comp.memberships(osmo$ppi, res)

# Compute (hard) clustering for nodes
res$clustering <- apply(res$comp.memb, 2, which.max)

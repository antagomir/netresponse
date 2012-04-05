# 2. netresponse test
# test later with varying parameters

# Load the package
library(netresponse)

data(toydata)
D <- toydata$emat
netw <- toydata$netw

# The toy data is random data with 10 features (genes). 
# The features 
rf <- c(4, 5, 6)
#form a subnetwork with coherent responses
# with means 
r1 <- c(0, 3, 0)
r2 <- c(-5, 0, 2)
r3 <- c(5, -3, -3)
mu.real <- rbind(r1, r2, r3)
# real weights
w.real <- c(70, 70, 60)/200
# and unit variances
rv <- 1

# Compute the model
#res <- detect.responses(D, netw, verbose = TRUE, mc.cores = 2)
res <- detect.responses(D, netw, verbose = TRUE, max.responses = 4)
# FIXME: we only get the correct 3-mode solution with
# max.responses <- 4; fix this.

# Subnets (each is a list of nodes)
subnet <- get.subnets(res)

# the correct subnet is retrieved in subnet number 2:
#> subnet[[2]]
#[1] "feat4" "feat5" "feat6"

# how about responses
# Retrieve model for the subnetwork with lowest cost function value
# means, standard devations and weights for the components
m <- get.model.parameters(res, subnet.id = "Subnet-2")  

# order retrieved and real response means by the first feature 
# (to ensure responses are listed in the same order)
# and compare deviation from correct solution
ord.obs <- order(m$mu[,1])
ord.real <- order(mu.real[,1])

print(paste("Correlation between real and observed responses:", cor(as.vector(m$mu[ord.obs,]), as.vector(mu.real[ord.real,]))))

# all real variances are 1, compare to observed ones
print(paste("Maximum deviation from real variances: ", max(abs(rv - range(m$sd))/rv)))

# weights deviate somewhat, this is likely due to relatively small sample size
#print("Maximum deviation from real weights: ")
#print( (w.real[ord.real] - m$w[ord.obs])/w.real[ord.real])

print("estimated and real mean matrices")
print(m$mu[ord.obs,])
print(mu.real[ord.real,])

# Test these later: some problems occur during build

# sample-response assignments for given subnet
# response.probabilities <- response2sample(res, subnet.id = 1)

# Internal. Should this be visible?
# means, standard devations and weights for the components
# for one subnet
# model <- get.model(res, subnet.id = 1) 

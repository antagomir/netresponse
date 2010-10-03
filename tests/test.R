# 2. netresponse test
# test later with varying parameters

# Load the package
library(netresponse)

data(toydata)
D <- toydata$emat
netw <- toydata$netw

# Compute the model
res <- detect.responses(D, netw, verbose = FALSE)

# Subnets (each is a list of nodes)
subnet <- get.subnets(res)


# Test these later: some problems occur during build

# sample-response assignments for given subnet
# response.probabilities <- response2sample(res, subnet.id = 1)

# Retrieve model for the subnetwork with lowest cost function value
# means, standard devations and weights for the components
#m <- get.model.parameters(res, subnet.id = 1)  

# Internal. Should this be visible?
# means, standard devations and weights for the components
# for one subnet
# model <- get.model(res, subnet.id = 1) 

# components for given subnet (nodes for each component)
# subnet.components(res, subnet.id = 1)
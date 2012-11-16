# Visualization

library(netresponse)

source("toydata.R")

######################################################################

res <- detect.responses(D, verbose = TRUE, max.responses = 10, mixture.method = "bic", information.criterion = "BIC", merging.threshold = 1, bic.threshold = 10)

# Subnets (each is a list of nodes)
subnets <- get.subnets(res)
subnet.id <- names(subnets)[[1]]

# plot subnetwork components
#vis <- plot.responses(res, subnet.id, plot.mode = "network")
vis <- plot.responses(res, subnet.id, plot.mode = "pca")
vis <- plot.responses(res, subnet.id, plot.mode = "heatmap")
vis <- plot.responses(res, subnet.id, plot.mode = "boxplot.data")
vis <- plot.responses(res, subnet.id, plot.mode = "response.barplot")


# plot subnetwork
#plot.subnet(res, subnet.id)

# plot color scale
#plot.scale(vis$breaks, vis$palette, two.sided = TRUE)

# Visualize data and centroids for one subnetwork with PCA
# pca.plot(res, subnet.id)


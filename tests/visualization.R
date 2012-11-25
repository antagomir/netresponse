# Visualization

library(netresponse)

source("toydata.R")

# --------------------------------------------------------------------

set.seed(4243)

res <- detect.responses(D, verbose = TRUE, max.responses = 10, mixture.method = "bic", information.criterion = "BIC", merging.threshold = 1, bic.threshold = 10, pca.basis = FALSE)

# Subnets (each is a list of nodes)
subnet.id <- names(get.subnets(res))[[1]]

# --------------------------------------------------------------------

# Plot subnetwork components
vis <- plot.responses(res, subnet.id, plot.mode = "pca")
vis <- plot.responses(res, subnet.id, plot.mode = "network")
vis <- plot.responses(res, subnet.id, plot.mode = "heatmap")
vis <- plot.responses(res, subnet.id, plot.mode = "boxplot.data")
vis <- plot.responses(res, subnet.id, plot.mode = "response.barplot")

# --------------------------------------------------------------------

# plot color scale
plot.scale(vis$breaks, vis$palette, two.sided = TRUE)



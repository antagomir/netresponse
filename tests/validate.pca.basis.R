
skip <- FALSE

if (!skip) {
# Visualization

library(netresponse)

#fs <- list.files("~/Rpackages/netresponse/netresponse/R/", full.names = T); for (f in fs) {source(f)}

source("toydata2.R")

# --------------------------------------------------------------------

set.seed(4243)
mixture.method <- "bic"

# --------------------------------------------------------------------

res <- detect.responses(D, verbose = TRUE, max.responses = 10, 
	   		       mixture.method = mixture.method, information.criterion = "BIC", 
			       merging.threshold = 1, bic.threshold = 10, pca.basis = FALSE)

res.pca <- detect.responses(D, verbose = TRUE, max.responses = 10, mixture.method = mixture.method, information.criterion = "BIC", merging.threshold = 1, bic.threshold = 10, pca.basis = TRUE)

# --------------------------------------------------------------------

k <- 1

# Incorrect VDP: two modes detected
# Correct BIC: single mode detected
subnet.id <- names(get.subnets(res))[[k]]

# Correct: single mode detected (VDP & BIC)
subnet.id.pca <- names(get.subnets(res.pca))[[k]]

# --------------------------------------------------------------------------------------------------

vis1 <- plot_responses(res, subnet.id, plot_mode = "pca", main = paste("NoPCA; NoDM"))
vis2 <- plot_responses(res, subnet.id, plot_mode = "pca", datamatrix = D, main = "NoPCA, DM")
vis3 <- plot_responses(res.pca, subnet.id.pca, plot_mode = "pca", main = "PCA, NoDM")
vis4 <- plot_responses(res.pca, subnet.id.pca, plot_mode = "pca", datamatrix = D, main = "PCA, DM")

# With original data: VDP overlearns; BIC works; with full covariance data 
# With PCA basis: modes detected ok with both VDP and BIC.

# ------------------------------------------------------------------------

# TODO
# pca.plot(res, subnet.id)
# plot_subnet(res, subnet.id) 
}
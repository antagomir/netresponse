#' PlotMixtureMultivariate.deprecated
#' 
#' Visualize data, centroids and response confidence intervals for a given
#' Gaussian mixture model with PCA. Optionally, color the samples according to annotations
#' labels.
#' 
#' @param x data matrix (samples x features)
#' @param means mode centroids (modes x features)
#' @param sds mode standard deviations, assuming diagonal covariance matrices (modes x features, each row giving the sqrt of covariance diagonal for the corresponding mode)
#' @param ws weight for each mode
#' @param labels Optional: sample class labels to be indicated in colors.
#' @param confidence Confidence interval for the responses based on the
#' covariances of each response. If NULL, no plotting.
#' @param ... Further arguments for plot function.
#' @return Used for its side-effects.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation('netresponse') for citation details.
#' @keywords utilities
#' @examples #plotMixture(dat, means, sds, ws)
PlotMixtureMultivariate.deprecated <- function(x, means, sds, ws, labels = NULL, 
    confidence = 0.95, ...) {
    
    if (!is.null(labels)) {
        if (is.null(names(labels))) {
            names(labels) <- rownames(x)
        }
    }
    
    # Add cluster centroids
    dat2 <- rbind(x, means)
    
    # center the data
    dat.mean <- colMeans(x)
    dat.centered <- t(t(dat2) - dat.mean)
    
    # PCA, two principal components
    pca <- princomp(dat.centered)
    
    # projection plane
    v <- as.matrix(pca$loadings[, 1:2])
    
    # Projected centroids (in PC space) for the detected components
    dat.pca <- dat.centered %*% v
    
    nlab <- length(unique(labels))
    if (nlab > 1) {
        my.palette <- palette(rainbow(nlab))
        cols <- my.palette[labels]
    } else {
        cols <- "black"
    }
    
    plot(dat.pca[1:nrow(x), ], main = paste("PCA plot"), xaxt = "n", yaxt = "n", 
        xlab = "", ylab = "", col = cols, pch = 19)
    
    if (nlab > 1) {
        legend("topleft", legend = as.character(unique(labels)), fill = unique(cols))
    }
    
    for (ri in 1:nrow(means)) {
        
        # Estimated covariance matrix for the response
        cmat <- diag(sds[ri, ]^2)
        
        # Projection of the covariance matrix in the PCA projection space
        cmat.projection <- diag(diag(t(v) %*% cmat %*% v))  # force it diagonal as it should be
        
        # Indicate estimated responses by ellipses
        if (!is.null(confidence)) {
            add.ellipse(centroid = dat.pca[ri, ], covmat = cmat.projection, confidence = confidence)
        }
    }
}


plot.associations <- function (x, subnet.id, labels, method = "hypergeometric", mode = "group.by.classes", ...) {
		  
  names(labels) <- rownames(x@datamatrix) # assumes that labels are in same order as data if names not given
  
  # Number of responses for this subnet
  N.responses <- length(x@models[[subnet.id]]$w)

  association.tab <- matrix(NA, nrow = length(unique(labels)), ncol = N.responses)
  rownames(association.tab) <- as.character(unique(labels))
  colnames(association.tab) <- names(x@models[[subnet.id]]$w)

  for (lab in unique(labels)) {
    # Pick samples associated with this annotation group
    samples <- names(labels)[which(labels == lab)]
    # Get enrichment information 
    # o <- order.responses(x, samples, method = "hypergeometric")$ordered.responses  # for all responses

    enr <- list()
    for (response in 1:N.responses) {
      enr[[response]] <- response.enrichment(subnet.id, x, samples, response, method)
    }

    neg.log.p <- -log10(sapply(enr, function (tab) {tab$info["pvalue"]}))
    association.tab[as.character(lab),] <- unname(neg.log.p)

  }

  if (mode == "group.by.responses") {
    barplot(association.tab, beside = TRUE, legend = TRUE, las = 1, cex.names = 0.8, main = "Label/response associations", ylab = "-log10(p)", ...)
  } else if (mode == "group.by.classes") {
    barplot(t(association.tab), beside = TRUE, legend = TRUE, las = 1, cex.names = 0.8, main = "Label/response associations", ylab = "-log10(p)", ...)
  }

}


plot.pca <- function (x, subnet.id, labels = NULL, confidence = 0.95, ...) {

  # FIXME: move these functions in NetResponse package

  if (!is.null(labels)) {
    if (is.null(names(labels))) {
       names(labels) <- rownames(x@datamatrix)
    }
  } 

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", "Subnet-", subnet.id, sep="")
  }
  
  # model parameters
  m <- get.model.parameters(x, subnet.id)
  
  # pick data (subnet only) 
  dat <- t(get.dat(x, subnet.id))

  # Add cluster centroids
  dat2 <- rbind(dat, m$mu)  

  # center the data
  dat.mean <- colMeans(dat)
  dat.centered <- t(t(dat2) - dat.mean)

  # PCA, two principal components
  pca <- princomp( dat.centered )

  # projection plane
  v <- as.matrix( pca$loadings[, 1:2] )
  
  # Projected centroids (in PC space) for the detected components
  dat.pca <- dat.centered %*% v

  nlab <- length(unique(labels))
  if (nlab > 1) {
    my.palette <- palette(rainbow(nlab)) 
    cols <- my.palette[labels]
  } else {
    cols <- "black"
  }

  plot(dat.pca[1:nrow(dat), ], main = paste("PCA plot: subnetwork ", subnet.id, sep = ""), xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", col = cols, pch = 19)

  if (nlab > 1) { legend("topleft", legend = as.character(unique(labels)), fill = unique(cols)) }

  for (rnam in rownames(m$mu)) { 

    # Estimated covariance matrix for the response
    cmat <- diag(m$sd[rnam,]^2)

    # Projection of the covariance matrix in the PCA projection space
    cmat.projection <- diag(diag(t(v)%*%cmat%*%v)) # force it diagonal as it should be

    # Indicate estimated responses by ellipses
    if (!is.null(confidence)){
      add.ellipse(centroid = dat.pca[rnam, ], covmat = cmat.projection, confidence = confidence)
    }
  }
}

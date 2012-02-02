# Copyright (C) 2010-2012 Leo Lahti
# Contact: Leo Lahti <leo.lahti@iki.fi>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.




# FIXME: make general plot.associations function
plotAssociations <- function (x, subnet.id, labels, method = "hypergeometric", mode = "group.by.classes", ...) {
		  
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


plotPCA <- function (x, subnet.id, labels = NULL, confidence = 0.95, ...) {

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

###########################################

plot.data <- function (x, subnet.id, labels, ...) {

    # ggplot2 boxplots for each user-defined sample category (listed in labels)
    library(ggplot2)
    dat <- t(netresponse::get.dat(x, subnet.id)) # samples x nodes
    df <- data.frame(list(labels = labels, dat))
    dfm <- melt(df, id = "labels")
    p <- ggplot(dfm) + aes(x = labels, y = dfm$value) + facet_wrap(~variable) + geom_boxplot() + opts(title = paste(subnet.id, ": annotation boxplot", sep = ""))    
    print(p)
    p

}

plot.expression <- function (x, maintext, ...) { # was: plot.matrix
		  
  # set color breakpoints and palette
  mybreaks <- set.breaks(1, interval = .02)
  mypalette <- colorRampPalette(c("blue", "black", "red"), space = "rgb")
  
  # compute differential expression in nodes with respect to the mean expression level for each gene
  ctrl.state <- colMeans(x)
  dmat <- t(t(x) - ctrl.state)

  # Color plot of the whole expression matrix, ordered by responses
  tmp <- plotMatrix.2way(dmat, mybreaks = mybreaks, maintext=maintext, cexlab=1, mypalette = mypalette)
}




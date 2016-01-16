#' @title Association strength between category labels and responses.
#' @description Plot association strength between user-defined category labels and responses in a selected subnetwork. Associations are showm in terms -log10(p) enrichment values for the annotation categories for the responses within the specified subnetwork. No correction for multiple testing. 
#' @param x NetResponseModel object
#' @param subnet.id Subnetwork.
#' @param labels Factor. Labels for the data samples. Name by samples, or
#' provide in the same order as in the original data.
#' @param method Method to calculate association strength.
#' @param mode group.by.responses or group.by.classes: indicate barplot
#' grouping type.
#' @param ... Other arguments to be passed for plot_
#' @return Used for side effect (plotting).
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso plot_responses
#' @references See citation("netresponse").
#' @keywords utilities
#' @examples #
plot_associations <- function (x, subnet.id, labels, 
		     	       method = "hypergeometric", 
			       mode = "group.by.classes", ...) {
		  
  # assumes that labels are in same order as data if names not given
  names(labels) <- rownames(x@datamatrix) 
  
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
 
      enr[[response]] <- response.enrichment(x[[subnet.id]]$qofz, 
      		      	 			samples, response, method)
    }    
    neg.log.p <- -log10(sapply(enr, function (tab) {tab$info["pvalue"]}))
    association.tab[as.character(lab),] <- unname(neg.log.p)

  }

  if (mode == "group.by.responses") {
    barplot(association.tab, beside = TRUE, legend = TRUE, las = 1, 
    	    cex.names = 0.8, main = "Label/response associations", 
	    ylab = "-log10(p)", ...)
  } else if (mode == "group.by.classes") {
    barplot(t(association.tab), beside = TRUE, legend = TRUE, las = 1, 
    	    cex.names = 0.8, main = "Label/response associations", 
	    ylab = "-log10(p)", ...)
  }

}




#' @title plotPCA
#' @description Visualize data, centroids and response confidence intervals for a given
#' subnetwork with PCA. Optionally, color the samples according to annotations
#' labels.
#' @param x NetResponseModel object. Output from the detect.responses function.
#' @param subnet.id Subnetwork id. Either character as 'Subnetwork-2' or
#' numeric as 2, which is then converted to character.
#' @param labels Optional: sample class labels to be indicated in colors.
#' @param confidence Confidence interval for the responses based on the
#' covariances of each response. If NULL, no plotting.
#' @param ... Further arguments for plot function.
#' @return Used for its side-effects.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse") for citation details.
#' @keywords utilities
#' @examples #plotPCA(x, subnet.id)
plotPCA <- function (x, subnet.id, labels = NULL, confidence = 0.95, ...) {

  if (!is.null(labels)) {
    if (is.null(names(labels))) {
       names(labels) <- rownames(x@datamatrix)
    }
  } 

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", 
    		       "Subnet-", subnet.id, sep="")
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

  plot(dat.pca[1:nrow(dat), ], main = paste("PCA plot: subnetwork ", 
  	subnet.id, sep = ""), xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", 
	col = cols, pch = 19)

  if (nlab > 1) { legend("topleft", 
     	     legend = as.character(unique(labels)), fill = unique(cols)) }

  for (ri in 1:nrow(m$mu)) { 

    # Estimated covariance matrix for the response
    cmat <- diag(m$sd[ri,]^2)

    # Projection of the covariance matrix in the PCA projection space
    # force it diagonal as it should be
    cmat.projection <- diag(diag(t(v)%*%cmat%*%v)) 

    # Indicate estimated responses by ellipses
    if (!is.null(confidence)){
      add.ellipse(centroid = dat.pca[ri, ], covmat = cmat.projection, 
         confidence = confidence)
    }
  }
}



#' @title PlotMixtureBivariate
#' @description Visualize data, centroids and response confidence intervals for a given
#' Gaussian mixture model in two-dimensional (bivariate) case. Optionally, 
#' color the samples according to annotations labels.
#' @param x data matrix (samples x features)
#' @param means mode centroids (modes x features)
#' @param sds mode standard deviations, assuming diagonal covariance 
#'    matrices (modes x features, each row giving the sqrt of covariance 
#'    diagonal for the corresponding mode)
#' @param ws weight for each mode
#' @param labels Optional: sample class labels to be indicated in colors.
#' @param confidence Confidence interval for the responses based on the
#' covariances of each response. If NULL, no plotting.
#' @param main title text
#' @param ... Further arguments for plot function.
#' @return Used for its side-effects.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse") for citation details.
#' @keywords utilities
#' @export
#' @examples #plotMixture(dat, means, sds, ws)
PlotMixtureBivariate <- function (x, means, sds, ws, labels = NULL, 
		     confidence = 0.95, main = "", ...) {

  if (!is.null(labels)) {
    if (is.null(names(labels))) {
       names(labels) <- rownames(x)
    }
  } 

  # Add cluster centroids
  dat.proj <- rbind(x, means)  
 
  nlab <- length(unique(labels))
  if (nlab > 1) {
    my.palette <- palette(rainbow(nlab)) 
    cols <- my.palette[labels]
  } else {
    cols <- "black"
  }

  plot(dat.proj[1:nrow(x), ], main = main, xaxt = 'n', yaxt = 'n', 
  			   xlab = "", ylab = "", col = cols, pch = 19)

  if (nlab > 1) { legend("topleft", legend = as.character(unique(labels)), 
     	     	  		    fill = unique(cols)) }

  for (ri in 1:nrow(means)) { 

    # Estimated covariance matrix for the response
    cmat <- diag(sds[ri,]^2)

    # Indicate estimated responses by ellipses
    if (!is.null(confidence)){
      add.ellipse(centroid = dat.proj[ri, ], covmat = cmat, 
      			   confidence = confidence)
    }
  }
}

#' @title PlotMixtureMultivariate
#' @description Visualize data, centroids and response confidence intervals for a given
#' Gaussian mixture model with PCA. Optionally, color the samples according 
#' to annotations labels.
#' @param x data matrix (samples x features)
#' @param means mode centroids (modes x features)
#' @param sds mode standard deviations, assuming diagonal covariance matrices 
#'        (modes x features, each row giving the sqrt of covariance diagonal 
#' 	  for the corresponding mode)
#' @param ws weight for each mode
#' @param labels Optional: sample class labels to be indicated in colors.
#' @param title title
#' @param modes Optional: provide sample modes for visualization already in 
#' 	  the input
#' @param pca The data is projected on PCA plane by default (pca = TRUE). 
#' 	  By setting this off (pca = FALSE) it is possible to visualize 
#' 	  two-dimensional data in the original domain.
#' @param qofz Sample-response probabilistic assignments matrix 
#' 	  (samples x responses)
#' @param ... Further arguments for plot function.
#' @return Used for its side-effects.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse") for citation details.
#' @keywords utilities
#' @export
#' @examples #plotMixture(dat, means, sds, ws)
PlotMixtureMultivariate <- function (x, means, sds, ws, labels = NULL, 
			title = NULL, modes = NULL, pca = FALSE, 
			qofz = NULL, ...) {

  # x <- t(X); means = model$params$mu; sds = model$params$sd; 
  # ws = model$params$w; labels <- NULL; title = ""

  # Circumvent warnings
  Comp.1 <- Comp.2 <- dat2 <- NULL

  if (!is.null(labels) && is.null(names(labels))) {
    names(labels) <- rownames(x)
  } 

  if (pca || ncol(x) > 2) {

      # center the data
      dat.mean <- colMeans(x)
      dat.centered <- t(t(x) - dat.mean)

      # PCA, two principal components
      pca <- princomp( dat.centered )

      # projection plane
      v <- as.matrix( pca$loadings[, 1:2] )
  
      # Projected centroids (in PC space) for the detected components
      dat.pca <- dat.centered %*% v

      pld <- dat.pca

      xtitle <- "PCA1"
      ytitle <- "PCA2"

      if (is.null(title)) {title <- paste("PCA (", ncol(x), " features)", 
      	 sep = "")}

  } else {

    pld <- x
    xtitle <- colnames(x)[[1]]
    ytitle <- colnames(x)[[2]]
    if (is.null(title)) { title <- "Cross-plot" } 

  }

  nlab <- length(unique(labels))
  if (nlab > 1) {
    my.palette <- palette(rainbow(nlab)) 
    cols <- my.palette[labels]
  } else {
    cols <- "black"
  }


  if (is.null(modes)) {
  
    # Determine the most likely cluster for each sample (-> hard clusters)
    if (is.null(qofz)) {
      qofz <- P.r.s(t(x), list(mu = means, sd = sds, w = ws), log = TRUE)
      rownames(qofz) <- rownames(x)
      colnames(qofz) <- paste("Mode", 1:ncol(qofz), sep = "-")
    }
    modes <- paste("Mode", apply(qofz, 1, which.max), sep = "-")

  }

  # Form data.fframe
  df <- data.frame(list(Comp.1 = pld[, 1], Comp.2 = pld[, 2]))
  df$mode <- factor(modes)

  theme_set(theme_bw(15))
  p <- ggplot(df, aes(x = Comp.1, y = Comp.2, colour = mode)) + geom_point() 
  p <- p + ggtitle(title) + xlab(xtitle) + ylab(ytitle)
  print(p)

  p

}

#' @title Plot observed data. 
#' @description Plotting tool for measurement data.
#' Produces boxplot for each feature in each annotation category for the
#' selected subnetwork.
#' @param x NetResponseModel object.
#' @param subnet.id Specify the subnetwork.
#' @param labels Annotation categories.
#' @param ... Further arguments for plot function.
#' @return ggplot2 plot object
#' @author Leo Lahti <leo.lahti@@iki.fi>
#' @seealso plot_responses
#' @references See citation("netresponse")
#' @keywords utilities
#' @export
#' @examples #
plot_data <- function (x, subnet.id, labels, ...) {

    # ggplot2 boxplots for each user-defined sample category (listed in labels)
    dat <- t(get.dat(x, subnet.id)) # samples x nodes
    df <- data.frame(list(labels = labels, dat))
    dfm <- melt(df, id = "labels")
    p <- ggplot(dfm) + aes(x = labels, y = dfm$value) + facet_wrap(~variable) 
    p <- p + geom_boxplot() 
    p <- p + ggtitle(paste(subnet.id, ": annotation boxplot", sep = ""))
    print(p)
    p

}



#' @title plot_expression
#' @description Plot expression matrix in color scale. For one-channel data; plot expression
#' of each gene relative to its mean expression level over all samples. Blue
#' indicates decreased expression and red indicates increased expression.
#' Brightness of the color indicates magnitude of the change. Black denotes no
#' change.
#' @param x samples x features matrix
#' @param maintext main title
#' @param ... optional arguments
#' @return Used for its side effects.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso \code{\link{plot_scale}}
#' @references See citation("netresponse").
#' @keywords utilities
#' @export
#' @examples #plot_expression(x)
plot_expression <- function (x, maintext, ...) { # was: plot_matrix
		  
  # set color breakpoints and palette
  mybreaks <- set.breaks(1, interval = .02)
  mypalette <- colorRampPalette(c("blue", "black", "red"), space = "rgb")
  
  # compute differential expression in nodes with respect to the mean 
  # expression level for each gene
  ctrl.state <- colMeans(x)
  dmat <- t(t(x) - ctrl.state)

  # Color plot of the whole expression matrix, ordered by responses
  tmp <- plot_matrix(dmat, mybreaks = mybreaks, maintext=maintext, cexlab=1, 
      	 		   mypalette = mypalette)
}


#' @title plot_subnet
#' @description Plot the given subnetwork.
#' @param x Result from NetResponse (detect.responses function).
#' @param subnet.id Subnet id.
#' @param network Original network used in the modelling.
#' @param plot_names Plot node names (TRUE) or indices (FALSE).
#' @param ... Further arguments for plot function.
#' @return Used for its side-effects. Returns a matrix that describes the
#' investigated subnetwork.
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao. Maintainer:
#' Leo Lahti <leo.lahti@@iki.fi>
#' @references L. Lahti et al.: Global modeling of transcriptional responses in
#' interaction networks. Submitted.
#' @keywords utilities
#' @export
#' @examples #
#' # res <- detect.responses(D, netw, verbose = FALSE)
#' # net <- plot_subnet(res, subnet.idx = 1)
plot_subnet <- function (x, subnet.id, network, plot_names = TRUE, ...) {

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", 
    		       "Subnet-", subnet.id, sep="")
  }

  subnet.nodes <- get.subnets(x)[[subnet.id]]
  mynet <- network[subnet.nodes, subnet.nodes]

  tmp <- plot_response(x = NULL, mynet, mybreaks = NULL, mypalette = NULL, 
      	 		 colors = FALSE, maintext = subnet.id)

  mynet

}

 
 #message("convert to matrix graph format")
 #myg <- new("graphAM", mynet, "undirected")
 #myg2 <-as(myg, "graphNEL") 




#' @title plot_response
#' @description Plot a specific transcriptional response for a given subnetwork.
#' TRUE, colors = TRUE, plot_type = "twopi", ...)
#' @param x A numerical vector, or NULL.
#' @param mynet Binary matrix specifying the interactions between nodes.
#' @param mybreaks Specify breakpoints for color plot_
#' @param mypalette Specify palette for color plot_
#' @param plot_names Plot node names (TRUE) or indices (FALSE).
#' @param colors Plot colors. Logical.
#' @param plot_type Network plot mode. For instance, 'neato' or 'twopi'.
#' @param ... Further arguments for plot function.
#' @return Used for its side-effects.
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao. Maintainer:
#' Leo Lahti <leo.lahti@@iki.fi>
#' @references L. Lahti et al.: Global modeling of transcriptional responses in
#' interaction networks. Submitted.
#' @keywords utilities
#' @export
#' @examples 
#' #tmp <- plot_response(model, mynet, 
#' #  	maintext = paste("Subnetwork", subnet.id))
plot_response <- function (x, mynet, mybreaks, mypalette, plot_names = TRUE, colors = TRUE, 
	     plot_type = "twopi", ...) {

  check.bins <- function (difexp, mybreaks) {

  # check color scale bin for each expression value
  bins <- c()
  for (i in 1:length(difexp)) {
      # which color bins are smaller than our difexp value
      # (for probet: i, mode:mode)
      inds <- which(difexp[[i]] > mybreaks)
      if (length(inds) == 0) {
        bins[[i]] <- 1
      } else if (length(inds) > 0)  {
        bins[[i]] <- max(inds) + 1
      }
  }

    bins
  }


   # Add node color for specific nodes
   nAttrs <- list()
   if (colors) {
     bins <- check.bins(x, mybreaks)
     nAttrs$fillcolor <- mypalette(length(mybreaks) + 1)[bins]
   } else {
     nAttrs$fillcolor <- rep("white", nrow(mynet))
   }
   names(nAttrs$fillcolor) <- rownames(mynet)

   # add node names for all nodes
   if (plot_names) {
     nodenames <- rownames(mynet)
   } else {
     nodenames <- rep("", nrow(mynet))
   }

   nAttrs$label <- nodenames
   names(nAttrs$label) <- rownames(mynet)

   myg <- as(new("graphAM", mynet, "undirected"), "graphNEL")
   
   plot(myg, y = plot_type, nodeAttrs = nAttrs, ...)

 }





#' @title plot_responses
#' @description Plot the detected transcriptional responses for a given subnetwork.
#' plot_mode = "network", xaxis = TRUE, yaxis = TRUE, plot_type = "twopi", mar
#' = c(5, 4, 4, 2), horiz = TRUE, datamatrix = NULL, scale = FALSE, ...)
#' @param x Result from NetResponse (detect.responses function).
#' @param subnet.id Subnet id.
#' @param nc Number of columns for an array of images.
#' @param plot_names Plot node names (TRUE) or indices (FALSE).
#' @param plot_mode network: plot responses as a subnetwork graph; matrix,
#' heatmap: plot subnetwork expression matrix. For both, expression of each
#' gene is shown relative to the mean expression level of the gene;
#' boxplot_data: feature-wise boxplots for hard sample-to-response assignments;
#' response.barplot: estimated response centroids as barplot including 95%
#' confidence intervals for the means; pca: PCA projection with estimated
#' centroids and 95% intervals. In 1-dimensional case a histogram is drawn. In
#' two-dimensional case the original coordinates are used.
#' @param xaxis,yaxis Logical. Plot row/column names.
#' @param plot_type Network plot mode. For instance, 'neato' or 'twopi'.
#' @param mar Figure margins.
#' @param horiz Logical. Horizontal barplot_
#' @param datamatrix datamatrix
#' @param scale scale the phylotypes to unit length (only implemented for plot_mode = "matrix"
#' @param ... Further arguments for plot function.
#' @return Used for its side-effects.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso \code{\link{plot_scale}}
#' @references See citation("netresponse")
#' @keywords utilities
#' @export
#' @examples #
#' #res <- detect.responses(D, netw)
#' #vis <- plot_responses(res, subnet.id)
plot_responses <- function (x, subnet.id, nc = 3, plot_names = TRUE, 
	                    plot_mode = "network", xaxis = TRUE, yaxis = TRUE, 
			    plot_type = "twopi", mar = c(5, 4, 4, 2), 
			    horiz = TRUE, datamatrix = NULL, 
			    scale = FALSE, ...) {

  # xaxis = TRUE; yaxis = TRUE; plot_type = "twopi"; mar = c(5, 4, 4, 2); 
  # horiz = TRUE; datamatrix = NULL; scale = FALSE; x <- res; nc <- 3; 
  # plot_names = TRUE; plot_mode = "pca"; main = paste("NoPCA; NoDM")

  responses <- NULL
  variable <- NULL
  p <- NULL # return ggplot2 object if available

  if (is.null(datamatrix)) {
    datamatrix <- x@datamatrix
  }

  value <- tmp <- NULL

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", "Subnet-", subnet.id, sep="")
  }

  pars <- get.model.parameters(x, subnet.id)
  subnet.nodes <- get.subnets(x)[[subnet.id]]
  datamatrix <- datamatrix[, subnet.nodes]

  if (class(x@network) %in% c("graphNEL", "igraph")) {
    # FIXME speedup by just converting the correct subset, not whole matrix
    mynet <- as(x@network, "graphAM")@adjMat
    rownames(mynet) <- colnames(mynet) <- nodes(x@network)
    mynet <- mynet[subnet.nodes, subnet.nodes]
  } else {
    mynet <- x@network[subnet.nodes, subnet.nodes]
  }

  # Set color breakpoints and palette
  mybreaks <- set.breaks(1, interval = .02)
  mypalette <- colorRampPalette(c("blue", "white", "red"), space = "rgb")
  
  # compute differential expression in nodes with respect to the mean 
  # expression level for each gene
  ctrl.state <- colMeans(datamatrix)
  centroids <- t(pars$mu)
  difexp <- apply(centroids, 2, function(x){ x - ctrl.state })
  rownames(difexp) <- rownames(centroids)

  # Determine the most likely mode for each sample (-> hard clusters)
  qofz <- x@models[[subnet.id]]$qofz
  modes <- apply(qofz, 1, which.max)

  if (plot_mode == "network") {
    par(mfrow = c(ceiling(length(pars$w)/nc), nc))
    for (comp in 1:length(pars$w)) {
      tmp <- plot_response(difexp[,comp], mynet, mybreaks, mypalette, 
      	     		   plot_names,
                           main = paste(subnet.id, "/Response-", comp, sep=""), 
			   plot_type = plot_type, ...)
    }
  } else if (plot_mode == "matrix" || plot_mode == "heatmap") {

    # order samples according to responses
    s2r <- apply(x[[subnet.id]]$qofz, 1, which.max) 
    ordered.samples <- order(s2r)

    dmat <- datamatrix[ordered.samples, subnet.nodes]
    dmat <- t(t(dmat) - ctrl.state)

    if (scale) {
      dmat <- scale(dmat, center = FALSE, scale = TRUE)
    }

    par(mfrow = c())

    # Color plot of the whole expression matrix, ordered by responses
    if (horiz) {
      tmp <- plot_matrix(t(dmat), mybreaks = mybreaks, maintext = subnet.id, xlab=NULL, ylab=NULL, palette = mypalette, xaxis = yaxis, yaxis = xaxis, mar = mar, ...)
    } else {
      tmp <- plot_matrix(dmat, mybreaks = mybreaks, maintext = subnet.id, 
      xlab=NULL, ylab=NULL, palette = mypalette, xaxis = xaxis, yaxis = yaxis, 
      mar = mar, ...)
    }

  } else if (plot_mode == "boxplot_data") {

    s2r <- apply(x[[subnet.id]]$qofz, 1, which.max)
    label <- factor(s2r)  

    # Ggplot2 boxplot handy as determines the grid size automatically
    # List samples in each response (hard assignments)
    dat <- t(get.dat(x, subnet.id)) # samples x nodes
    df <- data.frame(list(responses = label, dat))
    dfm <- melt(df, id = "responses")

    p <- ggplot(dfm) + aes(x = responses, y = value) + facet_wrap(~variable) 
    p <- p + geom_boxplot() 
    p <- p + ggtitle(paste(subnet.id, ": response boxplot", sep = ""))

    print(p)

  } else if (plot_mode == "response.barplot") {

    # FIXME: does not work

    # Plot cross-bars for estimated means and 95% intervals for each 
    # response for each node
    m <- get.model.parameters(x, subnet.id)

    # Use 1.96*std of the mean (95% quantile) for error limits
    # using the soft assignment sum as the sample size for each cluster
    # FIXME: use directly the parametric estimates from the model?
    # NOTE: that wouldn't work with pca.basis version directly

    # note: reverse the groups since also sds is reversed
    s2r <- apply(x[[subnet.id]]$qofz, 1, which.max)

    dat <- t(get.dat(x, subnet.id)) # samples x nodes
    response <- factor(s2r) 
    #factor(apply(sample2response(x, subnet.id), 1, which.max))
    df <- data.frame(dat)
    df$response <- response
    dfm <- melt(df, id.var = "response")    

    ggplot(dfm) + aes(x = response, y = value, fill = response) + 
    facet_wrap(~variable) + geom_bar(stat="identity")

    # mean and std of mean
    df <- ddply(dfm, c("response", "variable"), function (dd) {
       c(mean = mean(dd$value), 
         sd = 1.96*sd(dd$value)/sqrt(sum(x[[subnet.id]]$qofz[, dd$response])))})

    p <- qplot(response, mean, fill=variable, data=df, geom="bar", 
      	 		 position="dodge")
    p <- p + geom_errorbar(aes(ymax=mean+sd, ymin=mean-sd), 
      	     		 position="dodge") + theme_bw()

    print(p)


  } else if (plot_mode == "pca") {
    
    dmat <- datamatrix[, subnet.nodes]

    if (length(subnet.nodes) == 1) {

      p <- PlotMixtureUnivariate(dmat, 
      	    			 modes = modes,
				 title.text = subnet.id,
  				 xlab.text = "Component score", 
  				 ylab.text = "Frequency", 
  				 binwidth = 0.1)


    } else {

      if (ncol(dmat) > 2) {pca <- TRUE} else {pca <- FALSE}

      p <- PlotMixtureMultivariate(dmat, 
      	    			 modes = modes,
				 title = subnet.id, 
				 pca = pca)
					  
      # PlotMixtureBivariate(dmat, means, sds, ws, labels = NULL, 
      #				   confidence = 0.95, main = "") 
      # tmp <- plotPCA(x, subnet.id, labels = NULL, confidence = 0.95)

    }

  }

  list(breaks = mybreaks, palette = mypalette, info = tmp, p = p)

}



#' @title plot_scale
#' @description Plot the color scale used in visualization.
#' @param x Breakpoints for the plot_
#' @param y Color palette.
#' @param m Breakpoints' upper limit.
#' @param cex.axis Axis scale.
#' @param label.step Density of the labels.
#' @param interval Interval.
#' @param two.sided Plot two-sided (TRUE) or one-sided (FALSE) visualization.
#' @param label.start Label starting point.
#' @param Nlab Number of labels to plot_
#' @param ... Further arguments for plot function.
#' @return Used for its side-effects.
#' @author Leo Lahti <leo.lahti@@iki.fi>
#' @references See citation("netresponse")
#' @keywords utilities
#' @export
#' @examples #
#'   #res <- detect.responses(D, netw, verbose = FALSE)
#'   #vis <- plot_responses(res, subnet.idx)
#'   #plot_scale(vis$breaks, vis$palette)
plot_scale <- function (x, y, m = NULL, cex.axis = 1.5, label.step = 2, 
	      interval = .1, two.sided = TRUE, label.start = NULL, 
	      Nlab = 3, ...) {

  # x <- tmp$breaks; y <- tmp$palette; m = NULL; cex.axis = 1.5; 
  # label.step = 2; interval = .1; two.sided = TRUE; 
  # label.start = NULL; Nlab = 3

  if (two.sided) {
    
    if (length(m) > 0) {
      x <- set.breaks(m, interval)
    } else {
      mm <- max(x[-c(1, length(x))])
      m <- mm - interval/2
    }
  
    image(t(as.matrix(seq(-mm, mm, length = 100))), col = y(length(x) - 1),
          xaxt = 'n', yaxt = 'n', zlim = range(x), breaks = x)
    
    ndigits <- nchar(unlist(strsplit(as.character(mm), "\\."))[[2]])
    digit.step <- 10^(-ndigits)
    labs <- seq(-mm, mm, by = digit.step)
    label.start <- -max(abs(round(labs, ndigits)))
    start.position <- match(-label.start, round(labs, ndigits))
    end.position <- match(label.start, round(labs, ndigits))
    inds <- seq(start.position, end.position,length = Nlab)
      
    axis(2, at = inds/length(labs),
         labels = labs[inds], cex.axis = cex.axis, las = 2)
  }

  if (!two.sided) {

    mm <- max(x) + 1e6 # infty
    m <- max(x)
 
    labs <- seq(0, m, label.step)
    inds <- sapply(labs,function(lab){min(which(lab<=x))})
  
    image(t(as.matrix(seq(0, m, length = 100))), col = y(length(x) - 1),
          xaxt='n', yaxt='n', zlim = range(x), breaks = x)
    
    int <- 1/(length(x)-1)
    axis(2, at = seq(0, 1, by = int)[inds], labels = labs,
         cex.axis = cex.axis, las = 2)
  }
  
}



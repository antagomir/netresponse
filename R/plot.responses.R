plot.responses <-
function (x, subnet.id, nc = 3, plot.names = TRUE, plot.mode = "network", xaxis = TRUE, yaxis = TRUE, plot.type = "twopi", mar = c(5, 4, 4, 2), horiz = TRUE, ...) {

  require(igraph)
  require(Rgraphviz)

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", "Subnet-", subnet.id, sep="")
  }

  pars <- get.model.parameters(x, subnet.id)
  subnet.nodes <- get.subnets(x)[[subnet.id]]

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
  mypalette <- colorRampPalette(c("blue", "black", "red"), space = "rgb")
  
  # compute differential expression in nodes with respect to the mean expression level for each gene
  ctrl.state <- colMeans(x@datamatrix[, subnet.nodes])
  centroids <- t(pars$mu)
  difexp <- apply(centroids, 2, function(x){ x - ctrl.state })
  rownames(difexp) <- rownames(centroids)

  if (plot.mode == "network") {
    par(mfrow = c(ceiling(length(pars$w)/nc), nc))
    for (comp in 1:length(pars$w)) {
      tmp <- plot.response(difexp[,comp], mynet, mybreaks, mypalette, plot.names,
                          main = paste(subnet.id, "/Response-", comp, sep=""), plot.type = plot.type, ...)
    }
  } else if (plot.mode == "matrix" || plot.mode == "heatmap") {
    # order samples according to responses
    ordered.samples <- unlist(response2sample(x, subnet.id))

    dmat <- x@datamatrix[ordered.samples, subnet.nodes]
    dmat <- t(t(dmat) - ctrl.state)
    par(mfrow = c())

    # Color plot of the whole expression matrix, ordered by responses
    if (horiz) {
      tmp <- plotMatrix.2way(t(dmat), mybreaks = mybreaks, maintext = subnet.id, xlab="", ylab="", mypalette, xaxis = yaxis, yaxis = xaxis, mar = mar, ...)
    } else {
      tmp <- plotMatrix.2way(dmat, mybreaks = mybreaks, maintext = subnet.id, xlab="", ylab="", mypalette, xaxis = xaxis, yaxis = yaxis, mar = mar, ...)
    }


  } else if (plot.mode == "boxplot.data") {

    #r2s <- response2sample(x, subnet.id)
    label <- factor(apply(sample2response(x, subnet.id), 1, which.max))
    dat <- t(get.dat(x, subnet.id)) # samples x nodes
    nr <- ceiling(sqrt(ncol(dat)))
    nc <- ceiling(ncol(dat)/nr)
    par(mfrow = c(nr, nc))
    for (fnam in colnames(dat)) {
    	boxplot(dat[, fnam] ~ label, main = fnam)
    }
    tmp <- NULL
    
    # Ggplot2 boxplots, this is handy as it determines the grid size automatically
    # List samples in each response (hard assignments)
    #library(ggplot2)
    #s2r <- sample2response(x, subnet.id)
    #responses <- factor(apply(s2r, 1, which.max))
    #dat <- t(netresponse::get.dat(x, subnet.id)) # samples x nodes
    #df <- data.frame(list(responses = responses, dat))
    #dfm <- melt(df, id = "responses")
    #ggplot(dfm) + aes(x = responses, y = value) + facet_wrap(~variable) + geom_boxplot() + opts(title = paste(subnet.id, ": response boxplot", sep = ""))    


  } else if (plot.mode == "response.barplot") {

    # Plot cross-bars for estimated means and 95% intervals for each response for each node
    m <- get.model.parameters(x, subnet.id)

    # Flip the responses and features for visualization
    means <- m$mu[rev(1:nrow(m$mu)), rev(1:ncol(m$mu))]
    sds <- m$sd[rev(1:nrow(m$sd)), rev(1:ncol(m$sd))]

    par(mar = mar);
    bp <- seq(0.5, (nrow(means) + 1) * ncol(means), 1)
    # Use 1.96*std of the mean (95% quantile) for error limits
    # using the soft assignment sum as the sample size for each cluster
    # FIXME: later use directly the parametric estimates from the model?

    # note: reverse the groups since also sds is reversed
    std.of.mean <- sds/rev(sqrt(colSums(sample2response(x, subnet.id))))

    # Error bars
    eb1 <- rbind(rep(NA, ncol(means)), means - 1.96*std.of.mean)
    eb2 <- rbind(rep(NA, ncol(means)), means + 1.96*std.of.mean)

    # Plot error bars
    if (horiz) {
      tmp <- barplot(means, beside = TRUE, las = 1, cex.names = 0.8, legend = TRUE, main = subnet.id, horiz = horiz)
      segments(y0 = bp, x0 = eb1, y1 = bp, x1 = eb2, col = "gray20", lwd = 1.5)
    } else {
      tmp <- barplot(means, beside = TRUE, las = 2, cex.names = 0.8, legend = TRUE, main = subnet.id, horiz = horiz)
      segments(x0 = bp, y0 = eb1, x1 = bp, y1 = eb2, col = "gray20", lwd = 1.5)
    }

  }

  list(breaks = mybreaks, palette = mypalette, info = tmp)

}


#' filter.netw
#' 
#' Mostly for internal use. Prefilter edges if speedups required.
#' 
#' Include to the network only the edges with the highest mutual information,
#' calculated based on the first principal components.
#' 
#' @usage filter.netw(network, delta, datamatrix, params)
#' @param network network
#' @param delta associated cost function value changes for each node merge
#' @param datamatrix datamatrix
#' @param params parameters
#' @return Filtered network
#' @author Leo Lahti, Olli-Pekka Huovilainen and Antonio Gusmao.  Maintainer:
#' Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @examples #
filter.netw <- function (network, delta, datamatrix, params) {

  tmp <- list(network = network, delta = delta)

  if (params$max.subnet.size > 1) {
    if (params$verbose) {message("Filter the network to only keep the edges with highest mutual information")}
    # Filter out the least promising edges from the network
    # based on mutual information. For each variable, pick at most speedup.max.edges
    if (params$speedup && !is.null(params$speedup.max.edges)) {
      tmp <- filter.network(network, delta, datamatrix, params)
    } 
  }

  tmp

}



#' filter.network
#' 
#' 
#' Include to the network only the edges with the highest mutual information,
#' calculated based on the first principal components.
#' 
#' @param network network
#' @param delta associated cost function value changes for each node merge
#' @param datamatrix datamatrix
#' @param params parameters
#' @return Filtered network
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @examples #
filter.network <- function (network, delta, datamatrix, params) {

  # Include at maximum speedup.max.edges for each node in the network
  # based on mutual information
  # return list of edges to keep

  remove.edges <- c()
  uniq.edges <- unique(network[1,])  

  for (idx in 1:length(uniq.edges)) {
    if (params$verbose) {message(paste(idx, "/",  length(uniq.edges)))}
    a <- uniq.edges[[idx]]

    # Pick edges for this node
    eds <- which(network[1, ] == a)

    # Calculate MI scores
    #require(minet)
    mis <- c()
    mi.cnt <- 0  
    # FIXME: could be parallelized
    for (edge in eds){ # which(is.na(delta))
      mi.cnt <- mi.cnt + 1

      # Pick node indices
      i <- network[2, edge]

      # For singletons
      dat <- cbind(datamatrix[, a], datamatrix[, i])
                   
      mis[[mi.cnt]] <- build.mim(dat, estimator="mi.empirical", disc = "equalwidth", nbins = params$nbins)[1, 2]
    }
    
    # Edges to remove
    remove.edges <- eds[na.omit(order(mis, decreasing = TRUE)[-seq(params$speedup.max.edges)])]
    keep.edges <- setdiff(1:ncol(network), remove.edges)

    # Filter out remove.edges
    if (length(remove.edges) > 0) {
      network <- network[, keep.edges]
      delta <- delta[keep.edges] 
    }

  }

  list(network = network, delta = delta)

}


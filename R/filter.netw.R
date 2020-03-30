#' @title filter.netw
#' 
#' @description Mostly for internal use. Prefilter edges if speedups required.
#' 
#' @details Include only edges with the highest mutual information,
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
#' @references See citation('netresponse')
#' @keywords internal
filter.netw <- function(network, delta, datamatrix, params) {
    
    if (params$max.subnet.size > 1) {
        
        if (params$positive.edges) {
            tmp <- remove.negative.edges(network, delta, datamatrix)
            network <- tmp$network
            delta <- tmp$delta
        }
        
        # Filter out the least promising edges from the network based on mutual
        # information. For each variable, pick at most speedup.max.edges
        if (params$speedup && !is.null(params$speedup.max.edges)) {
            if (params$verbose) {
                message("Filter the network to only keep the edges with highest mutual information")
            }
            tmp <- filter.network(network, delta, datamatrix, params)
            network <- tmp$network
            delta <- tmp$delta
            
        }
        
    }
    
    list(network = network, delta = delta)
    
}




#' remove.negative.edges
#' 
#' Remove edges where the connected nodes correlate negatively
#' 
#' @param network network
#' @param delta associated cost function value changes for each node merge
#' @param datamatrix datamatrix
#' @return Filtered network
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation('netresponse')
#' @keywords internal
#' @examples #

remove.negative.edges <- function(network, delta, datamatrix) {
    
    message("Filtering out negative edges..")
    
    # Correlation matrix between the nodes FIXME: parallelize
    
    for (edge in seq_len(ncol(network))) {
        
        # Pick nodes for this edge
        a <- network[1, edge]
        b <- network[2, edge]
        
        # Calculate correlation btw. the two nodes
        cc <- cor(datamatrix[, a], datamatrix[, b], method = "spearman")
        
        # Give infinite cost for edges with non-positive association
        if (cc <= 0) {
            delta[[edge]] <- Inf
        }
    }
    
    # Remove edges with infinite cost
    keep.edges <- !is.infinite(delta)
    network <- network[, keep.edges]
    delta <- delta[keep.edges]
    
    list(network = network, delta = delta)
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
#' @references See citation('netresponse')
#' @keywords internal
#' @examples #
filter.network <- function(network, delta, datamatrix, params) {
    
    # Include at maximum speedup.max.edges for each node in the network based on
    # mutual information return list of edges to keep
    
    remove.edges <- c()
    uniq.edges <- unique(network[1, ])
    
    for (idx in seq_len(length(uniq.edges))) {
        if (params$verbose) {
            message(paste(idx, "/", length(uniq.edges)))
        }
        a <- uniq.edges[[idx]]
        
        # Pick edges for this node
        eds <- which(network[1, ] == a)
        
        # Calculate MI scores
        mis <- c()
        mi.cnt <- 0
        # FIXME: could be parallelized which(is.na(delta))
        for (edge in eds) {
            mi.cnt <- mi.cnt + 1
            
            # Pick node indices
            i <- network[2, edge]
            
            # For singletons
            dat <- cbind(datamatrix[, a], datamatrix[, i])
            
            mis[[mi.cnt]] <- build.mim(dat, estimator = "spearman")[1, 2]
        }
        
        # Edges to remove

        mis <- as.data.frame(lapply(mis, unlist))

        remove.edges <- eds[na.omit(order(mis, decreasing = TRUE)[-seq(params$speedup.max.edges)])]
        keep.edges <- setdiff(seq_len(ncol(network)), remove.edges)
        
        # Filter out remove.edges
        if (length(remove.edges) > 0) {
            network <- network[, keep.edges]
            delta <- delta[keep.edges]
        }
        
    }
    
    list(network = network, delta = delta)
    
}


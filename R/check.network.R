#' @title check.network
#' @description Internal use to check input network and format detect.responses.
#' @param network Input network, see detect.responses
#' @param datamatrix Input datamatrix, see detect.responses
#' @param verbose Print intermediate messages
#' @return 
#'   \item{formatted }{Formatted network (self-links removed)}
#'   \item{original }{Original network (possible in another representation format)} 
#'   \item{delta }{Cost function changes corresponding to the 'formatted' network.} 
#'   \item{nodes }{Nodes corresponding to the 'formatted' network.}
#' @author Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso detect.responses
#' @references See citation('netresponse')
#' @examples # check.network(network, datamatrix, verbose = FALSE)
check.network <- function(network, datamatrix, verbose = FALSE) {

    
    # If no network is given, assume fully connected net
    if (is.null(network)) {
        if (verbose) {
            warning("No network provided in function call: assuming fully connected nodes.")
        }
        network <- matrix(1, ncol(datamatrix), ncol(datamatrix))
        rownames(network) <- colnames(network) <- colnames(datamatrix)
    }

    # FIXME: later add other forms of sparse matrices from Matrix package
    accepted.formats.net <- c("matrix", "Matrix", "dgCMatrix", "dgeMatrix", "graphNEL", 
        "igraph", "graphAM")
    if (!is(network)[[1]] %in% accepted.formats.net) {
        stop(paste("network needs to be in one of the following formats:", paste(accepted.formats.net, 
            collapse = "; ")))
    }

    # Convert matrix into graphNEL if it is not already
    if (length(is(network)) > 1 || !("graphNEL" %in% is(network))) {

        if (any(is(network) %in% c("dgeMatrix", "dfCMatrix", "Matrix", "data.frame"))) {
            # Sparse matrix or data.frame needs to be converted first into matrix
            network <- as.matrix(network)
        }

        if (is.matrix(network)) {
            
            if (!nrow(network) == ncol(network)) {
                stop("Error: network nrow = ncol required.\n")
            }

            # Ensuring symmetric network
            if (any(!network == t(network))) {
                warning("Network is not symmetric. Removing link directions to force symmetric network.")
                network <- ((network + t(network)) > 0) - 0
            }
            
            # Filter non-connected nodes
            keep <- rowSums(network) > 1
            network <- network[keep, keep]
            
            # check that node names given in the data and correspond
            if (is.null(rownames(network)) || is.null(colnames(datamatrix))) {
                if (!nrow(network) == ncol(datamatrix)) {
                  stop("Error: Equal number of features required for the network and data matrix when feature names not given.\n")
                } else {
                  warning("Warning: network and/or data features are not named; matched by order.\n")
                  if (is.null(rownames(network)) && is.null(colnames(datamatrix))) {
                    rownames(network) <- colnames(network) <- as.character(seq_len(nrow(network)))
                    colnames(datamatrix) <- rownames(network)
                  } else if (is.null(rownames(network)) && !is.null(colnames(datamatrix))) {
                    rownames(network) <- colnames(network) <- colnames(datamatrix)
                  } else if (!is.null(rownames(network)) && is.null(colnames(datamatrix))) {
                    colnames(datamatrix) <- rownames(network)
                  }
                }
            }
            network <- as(new("graphAM", adjMat = network), "graphNEL")
        } else if (is(network) == "igraph") {
            network <- igraph.to.graphNEL(network)
        }
    }

    # Now the network is in graphNEL format
    
    # FIXME: adjust such that igraph does not need to be converted in graphNEL (which
    # is larger FIXME: add option to give this as input; seems to consume much less
    # memory than graphNEL

    # list network nodes that are not in datamatrix
    common.feats <- intersect(nodes(network), colnames(datamatrix))
    other.feats <- setdiff(nodes(network), common.feats)

    # remove features with no functional data
    if (length(other.feats) > 0) {
        if (verbose) {
            message(paste("removing network nodes that are not in datamatrix: ", 
                length(other.feats), " nodes removed;", length(common.feats), " nodes used for modeling."))
        }
        # message('converting to igraph')
        network <- igraph.from.graphNEL(network)
        # message('selecting nodes that have functional data')
        network <- subgraph(network, common.feats)
        # message('converting to graphNEL')
        network <- igraph.to.graphNEL(network)
        # network <- removeNode(other.feats, network)
    }
    
    message("convert the network into edge matrix")
    # store original network node list
    network.nodes <- nodes(network)
    network <- edgeMatrix(network, duplicates = FALSE)  # indices correspond to node list in network.nodes
    # order such that row1 < row2
    network <- apply(network, 2, sort)
    if (verbose) 
        message("removing self-links")
    network <- network[, !network[1, ] == network[2, ]]
    
    # Store the network in igraph format
    network.orig <- network  # store network used for modeling (preprocessed)
    # FIXME: igraph is more memory-efficient but could not be used as network is in
    # NetResponseModel definition for some reason. If possible, convert from graphNEL
    # to igraph later on.
    
    # Network rows: 1) nodes 2) nodes 3) merging cost function delta
    rownames(network) <- c("node1", "node2")
    
    # Delta corresponds to the columns of the 'network' object
    delta <- rep(NA, ncol(network))
    
    list(formatted = network, original = network.orig, delta = delta, nodes = network.nodes)
}

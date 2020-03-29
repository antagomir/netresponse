#' get.mis
#' 
#' Mainly for internal use. Estimate mutual information for node pairs based on
#' the first principal components
#' 
#' 
#' @usage get.mis(datamatrix, network, delta, network.nodes, G, params)
#' @param datamatrix datamatrix
#' @param network network
#' @param delta delta
#' @param network.nodes network.nodes
#' @param G G
#' @param params params
#' @return mutual information matrix
#' @author Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation('netresponse')
#' @keywords internal
#' @examples #
get.mis <- function(datamatrix, network, delta, network.nodes, G, params) {
    
    mis <- c()
    mi.cnt <- 0
    for (edge in which(is.na(delta))) {
        mi.cnt <- mi.cnt + 1
        # Pick node indices
        a <- network[1, edge]
        i <- network[2, edge]
        dat <- cbind(prcomp(matrix(datamatrix[, network.nodes[G[[a]]]], nrow(datamatrix)), 
            center = TRUE)$x[, 1], prcomp(matrix(datamatrix[, network.nodes[G[[i]]]], 
            nrow(datamatrix)), center = TRUE)$x[, 1])
        
        # mis[[mi.cnt]] <- build.mim(dat, estimator='mi.empirical', disc = 'equalwidth',
        # nbins = params$nbins)[1, 2]
        mis[[mi.cnt]] <- build.mim(dat, estimator = "spearman")[1, 2]
        
    }
    
    mis
}




build.mim <- function(dataset, estimator = "spearman", disc = "none", nbins = sqrt(NROW(dataset))) {
    # This function is licensed under cc-by-sa 3.0 Modified from minet 3.6.0 to use
    # discretize correctly (not exported in build.mim which may cause occasional
    # function name conflicts)
    
    # if (disc == 'equalfreq' || disc == 'equalwidth' || disc == 'globalequalwidth')
    # dataset <- infotheo::discretize(dataset, disc, nbins)
    
    if (estimator == "pearson" || estimator == "spearman" || estimator == "kendall") {
        mim <- cor(dataset, method = estimator, use = "complete.obs")^2
        diag(mim) <- 0
        maxi <- 0.999999
        mim[which(mim > maxi)] <- maxi
        mim <- -0.5 * log(1 - mim)
    } else {
        stop("unknown estimator")
    }
    # else if (estimator == 'mi.mm') estimator <- 'mm' else if (estimator ==
    # 'mi.empirical') estimator <- 'emp' else if (estimator == 'mi.sg') estimator <-
    # 'sg' else if (estimator == 'mi.shrink') estimator = 'shrink' if (estimator ==
    # 'mm' || estimator == 'emp' || estimator == 'sg' || estimator == 'shrink') { mim
    # <- mutinformation(dataset, method = estimator) diag(mim) <- 0 }
    mim[mim < 0] <- 0
    mim
}



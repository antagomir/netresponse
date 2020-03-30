#' @title Center data matrix.
#' @description Center data matrix to 0 for each variable by removing the means.
#' @param X The data set: samples x features. Each feature will be centered.
#' @param rm.na Ignore NAs.
#' @param meanvalue Can be used to set a desired center value. The default is 0.
#' @return Centered data matrix.
#' @note Note that the model assumes samples x features matrix, and centers
#' each feature.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation('netresponse')
#' @keywords utilities maths
#' @export
centerData <- function(X, rm.na = TRUE, meanvalue = NULL) {
    
    # Shift data matrix (columns) to zero, or given 'meanvalue'
    
    if (!rm.na) {
        xcenter <- colMeans(X)
        X2 <- X - rep(xcenter, rep.int(nrow(X), ncol(X)))
    } else {
        X2 <- array(NA, dim = c(nrow(X), ncol(X)), dimnames = dimnames(X))
        for (i in seq_len(ncol(X))) {
            x <- X[, i]
            nainds <- is.na(x)
            xmean <- mean(x[!nainds])
            X2[!nainds, i] <- x[!nainds] - xmean
        }
        dimnames(X2) <- dimnames(X)
    }
    
    if (!is.null(meanvalue)) {
        # Shift the data so that mean gets a specified value
        X2 <- X2 + meanvalue
    }
        
    X2
    
}

# The best academic advice I ever got was: 'Spend at least an hour every day on
# the manuscript closest to publication'.

#' @title generate.toydata
#' @description Generate toy data.
#' @usage D <- generate.toydata()
#' @param Dim Dimensionality of data
#' @param Nc Number of modes
#' @param Ns Number of data points
#' @param sd0 Component spread
#' @param rgam.shape Shape parameter for Gamma distribution 
#' @param rgam.scale Scale parameter for Gamma distribution
#' @return Simulated data matrix (samples x features)
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation('netresponse')
#' @keywords internal
#' @examples D <- generate.toydata()
#' @export
generate.toydata <- function(Dim = 3, Nc = 3, Ns = 200, sd0 = 3, rgam.shape = 2, 
    rgam.scale = 2) {
    
    # Generate means and variances (covariance diagonals) for the components
    component.means <- matrix(rnorm(Nc * Dim, mean = 0, sd = sd0), nrow = Nc, ncol = Dim)
    component.vars <- matrix(1/rgamma(Nc * Dim, shape = rgam.shape, scale = rgam.scale), 
        nrow = Nc, ncol = Dim)
    component.sds <- sqrt(component.vars)
    
    # Size for each component -> sample randomly for each data point from uniform
    # distr.  i.e. cluster assignments
    sample2comp <- sample.int(Nc, Ns, replace = TRUE)
    D <- array(NA, dim = c(Ns, Dim))
    for (i in seq_len(Ns)) {
        # component identity of this sample
        ci <- sample2comp[[i]]
        cm <- component.means[ci, ]
        csd <- component.sds[ci, ]
        D[i, ] <- rnorm(Dim, mean = cm, sd = csd)
    }
    
    colnames(D) <- paste("Feature-", seq_len(ncol(D)), sep = "")
    rownames(D) <- paste("Sample-", seq_len(nrow(D)), sep = "")
    
    list(data = D, means = component.means, sds = component.sds, sample2comp = sample2comp)
}

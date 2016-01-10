#' @title Univariate BIC mixture
#' @description Latent class analysis based on (infinite) Gaussian mixture
#' model. If the input (dat) is data matrix, a multivariate model is fitted. If
#' the input is a vector or a 1-dimensional matrix, a univariate model is
#' fitted.
#' @param x  dat vector (for univariate analysis) or a matrix (for multivariate analysis)
#' @param max.modes Maximum number of modes to be checked for mixture model selection
#' @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture.
#' @param min.modes minimum number of modes
#' @param ... Further optional arguments to be passed
#' @return Fitted latent class model (parameters and free energy)
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
bic.mixture.univariate <- function (x, max.modes, bic.threshold = 0, min.modes = 1, ...) { 

  # x <- datamatrix[, node];  max.modes = params$max.responses; bic.threshold = params$bic.threshold

  best.mode <- bic.select.best.mode(x, max.modes, bic.threshold, min.modes = min.modes) 
  mcl <- Mclust(x, G = best.mode)

  means <- as.vector(mcl$parameters$mean)
  sds <- as.vector(sqrt(mcl$parameters$variance$sigmasq))
  if (length(sds) == 1) {sds <- rep(sds, length(means))} 
  ws <- as.vector(mcl$parameters$pro)

  if (is.null(ws)) {warning("NULL weights, replacing with 1"); ws <- 1} 
  if (is.null(means)) {warning("NULL means, replacing with 1"); means <- 1} 
  if (is.null(sds)) {warning("NULL sds, replacing with 1"); sds <- 1} 

  Nparams <- length(means) + length(sds) + length(ws) 

  means <- matrix(means, nrow = length(ws))
  sds <- matrix(sds, nrow = length(ws))

  # Determine the most likely mode for each sample (-> hard clusters)
  # save(means, sds, ws, x, file = "~/tmp/tmp.RData")
  qofz <- P.r.s(matrix(x, nrow = 1), list(mu = means, sd = sds, w = ws), log = FALSE)
  rownames(qofz) <- names(x)

  names(means) <- names(sds) <- names(ws) <- paste("Mode", 1:length(ws), sep = "-")

  list(means = means, sds = sds, ws = ws, Nparams = Nparams, free.energy = -mcl$loglik, qofz = qofz)

}


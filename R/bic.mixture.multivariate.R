#' @title Multivariate BIC mixture
#' @description Latent class analysis based on (infinite) Gaussian mixture model. If the input (dat) is data matrix, a multivariate model is fitted. 
#' @param x matrix (for multivariate analysis)
#' @param max.modes Maximum number of modes to be checked for mixture model selection
#' @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture.
#' @param min.modes Minimum number of modes to be checked for mixture model selection
#' @param ... Further optional arguments to be passed
#' @return Fitted latent class model (parameters and free energy)
#' @export
#' @references See citation("netresponse") 
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
bic.mixture.multivariate <- function (x, max.modes, bic.threshold = 0, min.modes = 1, ...) { 

  # x <- mat; max.modes = params$max.responses; bic.threshold = params$bic.threshold

  best.mode <- bic.select.best.mode(x, max.modes, bic.threshold, min.modes) 

  mcl <- Mclust(x, G = best.mode)

  bic <- try(-mclustBIC(x, G = best.mode)[, "VVV"]) 
  if ( is.na(bic) ) { bic <- Inf } # infinitely bad = Inf

  means <- t(mcl$parameters$mean)
  vars <- t(apply(mcl$parameters$variance$sigma, 3, function(x){diag(x)}))
  sds <- sqrt(vars)
  ws <- as.vector(mcl$parameters$pro)
  if (is.null(ws)) {ws <- 1} 

  Nparams <- prod(dim(means)) + prod(dim(sds)) + length(ws) 

  # Determine the most likely mode for each sample (-> hard clusters)
  qofz <- P.r.s(t(x), list(mu = means, sd = sds, w = ws), log = FALSE)
  rownames(qofz) <- rownames(x)
  colnames(qofz) <- paste("Mode", 1:ncol(qofz), sep = "-")

  rownames(means) <- rownames(sds) <- names(ws) <- paste("Mode", 1:length(ws), sep = "-")
  colnames(means) <- colnames(sds) <- colnames(x)

  list(means = means, sds = sds, ws = ws, Nparams = Nparams, free.energy = -mcl$loglik, qofz = qofz, bic = bic)

}

#' @title Select best mode with BIC
#' @description Select optimal number of mixture components by adding components until 
#' the increase in objective function is below threshold.
#' @param x  dat vector (for univariate analysis) or a matrix (for multivariate analysis)
#' @param max.modes Maximum number of modes to be checked for mixture model selection
#' @param bic.threshold BIC threshold which needs to be exceeded before a new mode is added to the mixture.
#' @param min.modes Optiomal. Minimum number of modes.
#' @return Fitted latent class model (parameters and free energy)
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
bic.select.best.mode <- function (x, max.modes, bic.threshold, min.modes = 1) {

  # Cost for single mode
  # BIC : smaller is better
  # mclustBIC returns the value for -BIC, to be exact
  nc <- min.modes
  if (is.vector(x)) { # univariate
    m <- -mclustBIC(x, G = nc)[, "V"] 
  } else { # multivariate
    m <- -mclustBIC(x, G = nc)[, "VVV"] # BIC : smaller is better
  }
  
  # ----------------------------------------------------------------
  
  add.component <- TRUE
  best.mode <- min.modes
  if (max.modes == min.modes) {
    add.component <- FALSE
  }

  while (add.component && nc < max.modes) {

    nc <- nc + 1

    # BIC : smaller is better
    if (is.vector(x)) { # univariate
      m.new <- try(-mclustBIC(x, G = nc)[, "V"]) 
    } else { # multivariate
      m.new <- try(-mclustBIC(x, G = nc)[, "VVV"]) 
    }
    if ( is.na(m.new) ) { m.new <- Inf } # infinitely bad = Inf

    # FIXME: compressing data with PCA after dimensionality gets otherwise too high?
    # with around ncol(x) = 30 the mclustBIC is starting to produce NAs

    # FIXME: remove this when code works ok
    # if (is.na(m.new)) {save(x, nc, file = "m.new.RData")}
    
    bic.delta <- m.new - m

    if (bic.delta < -bic.threshold) { 
      best.mode <- nc 
      m <- m.new
    } else {
      add.component <- FALSE
    }
  }

  best.mode

}



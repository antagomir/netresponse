#' @title Plot univariate mixtures
#' @description Visualize data, centroids and stds for a given univariate Gaussian mixture model with PCA.
#' @param x data vector
#' @param means mode centroids
#' @param sds mode standard deviations
#' @param ws weight for each mode
#' @param title.text Plot title
#' @param xlab.text xlab.text
#' @param ylab.text ylab.text
#' @param binwidth binwidth for histogram
#' @param qofz Mode assignment probabilities for each sample. Samples x modes.
#' @param density.color Color for density lines
#' @param cluster.assignments Vector of cluster indices, indicating cluster for each data point
#' @param ... Further arguments for plot function.
#' @return Used for its side-effects
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse") for citation details.
#' @keywords utilities
#' @examples # plotMixtureUnivariate(dat, means, sds, ws)
PlotMixtureUnivariate <- function (x, means = NULL, sds = NULL, ws = NULL, title.text = NULL, xlab.text = NULL, ylab.text = NULL, binwidth = 0.05, qofz = NULL, density.color = "darkgray", cluster.assignments = NULL, ...) {
				 
  # Circumvent warnings
  ..density.. <- NULL
  vals <- NULL
  varname <- NULL
		 
  # Find cluster for each sample
  if (is.null(cluster.assignments)) {
    if (is.null(qofz)) {
      qofz <- P.r.s(t(matrix(x)), list(mu = means, sd = sds, w = ws), log = TRUE)
    }
    cluster.assignments <- apply(qofz, 1, which.max)
  }

  x <- unname(x)
  df <- data.frame(list(x = x))
  df$mode <- factor(cluster.assignments)

  # Histogram and density plot
  pg <- ggplot(df, aes(x=x)) 
  pg <- pg + geom_histogram(aes(fill = mode), binwidth=binwidth) 
  pg <- pg + theme_bw() + xlab(xlab.text) + ylab(ylab.text) 
  pg <- pg + ggtitle(title.text)

  # If normal mixture parameters are given, overlay Gaussians on top of the plot
  if (!is.null(means)) {

    # Estimated normal distributions from the mixture model
    h <- hist(x, seq(min(x) - binwidth - 1/binwidth, max(x) + binwidth + 1/binwidth, binwidth), plot = FALSE)
    df$vals <- seq(min(df$x), max(df$x), length=nrow(df)) # estimation points for fitted Gaussians
    scal <- max(h$counts)/max(ws[[1]]*dnorm(df$vals, mean = means[[1]], sd = sds[[1]]))
    for(comp in 1:length(means)){
      df2 <- data.frame(list(vals = df$vals, varname = scal*ws[[comp]]*dnorm(df$vals, mean = means[[comp]], sd = sds[[comp]])))
      # Add in normal distribution
      pg <- pg + geom_line(aes(x=vals, y=varname), color=density.color, data = df2) 
    }
  }

  pg

}



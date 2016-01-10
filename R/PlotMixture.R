#' @title Plot mixtures
#' @description Plot mixtures.
#' @param x data vector
#' @param qofz Mode assignment probabilities for each sample. Samples x modes.
#' @param binwidth binwidth for histogram
#' @param xlab.text xlab.text
#' @param ylab.text ylab.text
#' @param title.text title.text
#' @return Used for its side-effects
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse") for citation details.
#' @keywords utilities
#' @examples # PlotMixture(x, qofz)
PlotMixture <- function (x, qofz, binwidth = 0.05, xlab.text = NULL, ylab.text = NULL, title.text = NULL) {

  # hard sample clusters 		 
  sms <- apply(qofz, 1, which.max)
     
  x <- unname(x)
  df <- data.frame(list(x = x))
  df$mode <- factor(sms)

  # Histogram and density plot
  pg <- ggplot(df, aes(x=x)) 
  pg <- pg + geom_histogram(aes(fill = mode), binwidth = binwidth) 
  pg <- pg + theme_bw() + xlab(xlab.text) + ylab(ylab.text) 
  pg <- pg + ggtitle(title.text)

  pg

}

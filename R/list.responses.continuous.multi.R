#' @title Investigate association of a continuous variable and the modes
#' @description Investigate association of a continuous variable and the modes given a list of groupings
#' @param annotation.df annotation data.frame with discrete factor levels, rows
#' named by the samples
#' @param groupings Sample mode information. Each element corresponds to one grouping; each grouping lists samples for the modes within that grouping.
#' @param method method for quantifying the association
#' @param pth p-value threshold applied to adjusted p-values
#' @param verbose verbose 
#' @param rounding rounding digits
#' @return Table listing all associations between the factor levels and responses
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @importFrom qvalue qvalue
#' @keywords utilities
list.responses.continuous.multi <- function (annotation.df, groupings, method = "t-test", pth = Inf, verbose = TRUE, rounding = NULL) {

  #annotation.df <- annot[, continuous.vars]; groupings <- groupings.listing; pth = pth; method <- "t-test"; verbose = TRUE; rounding = NULL

  tab <- NULL				
  for (gn in names(groupings)) {

    gntab <- list.responses.continuous.single(annotation.df, groupings[[gn]], method = method, pth = Inf, verbose = verbose, rounding = rounding, adjust.p = FALSE) 

    gntab <- cbind(grouping = gn, gntab)
    tab <- rbind(tab, gntab)

  }

  tab$p.adj <- rep(NA, nrow(tab))
  if (nrow(tab) > 100) {
    qv <- qvalue(tab$pvalue, pi0.method = "bootstrap", fdr.level = 0.25)
    if (("qvalues" %in% names(qv))) {
      tab$p.adj <- qv$qvalues
    }
  } else {
    tab$p.adj <- p.adjust(tab$pvalue, method = "BH")
  }
  
  tab$pvalue <- NULL
  tab <- tab[tab$p.adj < pth,]
 
  tab

}

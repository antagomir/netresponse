#' @title List factor responses (minimal)
#' @description List significantly associated responses for all factors and levels in the given annotation matrix
#' @param annotation.df annotation data.frame with discrete factor levels, rows
#'   named by the samples
#' @param groupings List of groupings. Each model should have a sample-cluster assignment matrix qofz, or a vector of cluster indices named by the samples.
#' @param method method for enrichment calculation
#' @param min.size minimum sample size for a response
#' @param pth p-value threshold; applied to adjusted p-value
#' @param verbose verbose 
#' @param data data (samples x features; or a vector in univariate case) 
#' @param rounding rounding digits 
#' @return A list with two elements: Table listing all associations between the factor levels and responses; multiple p-value adjustment method
#' @author Contact: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @export
#' @keywords utilities
list.responses.factor.minimal <- function (annotation.df, groupings, method = "hypergeometric", min.size = 2, pth = Inf, verbose = TRUE, data = NULL, rounding = NULL) {

  # samples x features
  if(is.vector(data)) {
    data2 <- matrix(data)
    rownames(data2) <- names(data)
    data <- data2
  }

  if (!is.data.frame(annotation.df)) {
    stop("Provide data.frame for the annotation.df argument!")
  }

  # Collect the tables from all factors and levels here
  collected.table <- NULL

  for (fnam in colnames(annotation.df)) { 

    if (verbose) { message(fnam) }

    annotation.vector <- annotation.df[[fnam]]
    names(annotation.vector) <- rownames(annotation.df)

    # Order responses for each level of this factor
    responses.per.level <- factor.responses.minimal(annotation.vector, groupings, method = method, min.size = min.size, data = data) 
    responses.per.level <- responses.per.level[na.omit(names(responses.per.level))]

    # Add factor/level information
    mat <- NULL
    if (length(responses.per.level) > 0 && !is.null(na.omit(names(responses.per.level)))) {
      for (level in na.omit(names(responses.per.level))) {

        responses.per.level[[level]] <- cbind(responses.per.level[[level]], fnam, level)
        tmp <- responses.per.level[[level]]
        colnames(responses.per.level[[level]]) <- c(colnames(tmp)[1:(ncol(tmp)-2)], "Factor", "Level")

      }
    
      # Combine level-wise matrices
      mat <- responses.per.level[[1]]; 
      if (length(responses.per.level) > 1) {
        for (i in 2:length(responses.per.level)) { mat <- rbind(mat, responses.per.level[[i]])} 
      }

    }

    # FIXME: do.call(rbind, l) tai do.call(cbind, l)
    # is a neater way to implement this
    collected.table <- rbind(collected.table, mat)

  }

  # Order by pvalue
  collected.table <- collected.table[order(as.numeric(as.character(collected.table[, "pvalue"]))), ]
  collected.table <- data.frame(collected.table)
  collected.table$pvalue <- as.numeric(as.character(collected.table$pvalue))

  if (!is.null(collected.table)) {
  
    if (nrow(collected.table)>100) {
      collected.table$p.adj <- qvalue(collected.table$pvalue, gui = FALSE, fdr.level = 0.25)$qvalue
      p.adj.method <- "qvalue"
    } else {
      collected.table$p.adj <- p.adjust(collected.table$pvalue, method = "BH")
      p.adj.method <- "BH"
    }  

    # Significance Filtering 
    if (!all(is.na(collected.table$p.adj))) {
      collected.table <- collected.table[collected.table$p.adj < pth, ]
    } else {
      collected.table <- collected.table[collected.table$pvalue < pth, ]
    }
  }

  collected.table$mode <- as.character(collected.table$mode)
  collected.table$annotated.in.subset <- as.numeric(as.character(collected.table$annotated.in.subset))
  collected.table$fraction.in.subset <- as.numeric(as.character(collected.table$fraction.in.subset))
  collected.table$fraction.in.data <- as.numeric(as.character(collected.table$fraction.in.data))
  collected.table$pvalue <- NULL
  collected.table$qvalue <- NULL
  collected.table$p.adj <- as.numeric(as.character(collected.table$p.adj))
 
  if (!is.null(rounding)) {

    collected.table$fraction.in.subset <- round(collected.table$fraction.in.subset, rounding)
    collected.table$fraction.in.data <- round(collected.table$fraction.in.data, rounding)
    collected.table$p.adj <- round(collected.table$p.adj, rounding)
 
  }

  list(table = collected.table, p.adj.method = p.adj.method)

}


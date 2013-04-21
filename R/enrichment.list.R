enrichment.list <- function (groupings, annotation.vector) {

  annotated.samples <- names(which(!is.na(annotation.vector)))
  annotation.data <- annotation.vector[annotated.samples]
  names(annotation.data) <- annotated.samples

  associations <- NULL  	 
  for (sn in names(groupings)) {

    mode.samples <- groupings[[sn]]

    # annotated samples in the mode
    s <- intersect(mode.samples, annotated.samples)

    # annotated samples in other modes
    sc <- intersect(unlist(groupings[setdiff(names(groupings), sn)]), annotated.samples)

    if (length(na.omit(s)) > 1 && length(na.omit(sc)) > 1) {      
        pvals <- t.test(annotation.data[s], annotation.data[sc])$p.value
        fold.change <- mean(annotation.data[s]) - mean(annotation.data[sc])
     } else {
         warning(paste("Not enough annotated observations for response", sn))
         pvals <- NA
	 fold.change <- NA
     }

    associations <- rbind(associations, cbind(mode = sn, pvalue = pvals, fold.change = fold.change))
 
  }

  associations <- data.frame(list(mode = as.character(associations[, "mode"]), 
    		    	      pvalue = as.numeric(associations[, "pvalue"]), 
    		fold.change = as.numeric(associations[, "fold.change"])
    				    ), stringsAsFactors = FALSE)


  associations

}

enrichment.list <- function (models, annotation.vector) {

  annotated.samples <- names(which(!is.na(annotation.vector)))
  annotation.data <- annotation.vector[annotated.samples]
  names(annotation.data) <- annotated.samples

  if (is.null(names(models))) {names(models) <- 1:length(models)}

  associations <- NULL  	 
  for (sn in names(models)) {

    # samples in each mode (hard assignment)
    r2s <- response2sample(models[[sn]])

    pvals <- c()
    fold.change <- c()
    for (mo in 1:length(r2s)) {

      # annotated samples in the mode
      s <- intersect(r2s[[mo]], annotated.samples)

      # annotated samples in other modes
      sc <- intersect(unlist(r2s[-mo]), annotated.samples)

      if (length(na.omit(s)) > 1 && length(na.omit(sc)) > 1) {      
        pvals[[mo]] <- t.test(annotation.data[s], annotation.data[sc])$p.value
        fold.change[[mo]] <- mean(annotation.data[s]) - mean(annotation.data[sc])

       } else {
         warning(paste("Not enough annotated observations for response", mo))
         pvals[[mo]] <- NA
	 fold.change[[mo]] <- NA
       }
     }

     associations <- rbind(associations, cbind(subnet = rep(sn, length(r2s)), mode = 1:length(r2s), pvalue = pvals, fold.change = fold.change))

    }

  associations <- data.frame(list(model = associations[, "subnet"], 
    		    	      mode = associations[, "mode"], 
    		    	      pvalue = as.numeric(associations[, "pvalue"]), 
    		fold.change = as.numeric(associations[, "fold.change"])
    				    ))


  associations

}

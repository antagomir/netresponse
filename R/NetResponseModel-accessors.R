setMethod("get.subnets", "NetResponseModel", function (model, get.names = TRUE, min.size = 2, max.size = Inf, min.responses = 2) {

  #  Copyright (C) 2008-2011 Leo Lahti
  #  Licence: GPL >=2
                
  grouping <- model@last.grouping
  
  # Use feature names instead of indices
  if ( get.names ) {
    grouping <- lapply(grouping, function(x) {model@nodes[unlist(x)]})
  }
  
  # name the subnetworks
  names(grouping) <- paste("Subnet", 1:length(grouping), sep = "-")
        
  # If filters are given, apply them (stat needs to be specified)
	    
  # SUBNET SIZE
	        
  subnet.size <- sapply(grouping, length)    
  df <- data.frame(sapply(grouping, length))
  colnames(df) <- c("subnet.size")
  inds <- rownames(subset(df,
		   subnet.size >= min.size &
                   subnet.size <= max.size))
			    
  ## NUMBER OF RESPONSES
			        
  if ( min.responses > 1 ) {
				  
    stat <- model.stats( model )
        
    # check which subnets pass the filter
    inds.size <- rownames(subset(stat,
		          subnet.size >= min.size & subnet.size <= max.size))# & 
   
    inds.nresp <- rownames(stat)[which(stat[["subnet.responses"]] >= min.responses)]

    inds <- intersect(inds.size, inds.nresp)
    
  }
        
  # Get the filtered subnetwork list
  if (length(inds) == 0) {grouping <- NULL} else { grouping <- grouping[inds] }

  grouping
  	  
})


setMethod(f = "[[", signature("NetResponseModel"),
   definition = (function(x, i, j = "missing", ..., exact = TRUE) {
      if ( typeof(i) == "numeric" ){ i <- names(x)[[i]] }
      get.model.parameters(x, subnet.id = i)
   })
)


      #new("rpa.list", list(d = x$d[i,], sigma2 = x$sigma2[[i]], cind = x$cind, set = x$sets[[i]]))
#setReplaceMethod(f="[[",signature("ChromosomeArmModels"),
#                                definition=(function(x,i,j,value) {
#                                        x@models[[i]] <- value
#                                        return(x)
#                                }
#))
 

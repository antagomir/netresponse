get.subnets <- function (model, level = NULL, get.names = TRUE, stat = NULL, min.size = NULL, max.size = NULL, min.responses = NULL) {

  #  Copyright (C) 2008-2010 Leo Lahti
  #  Licence: GPL >=2
            
  # If no level specified, show the last iteration step
  if ( is.null(level) ) {
    grouping <- model@last.grouping
  } else { NULL }
  # FIXME: add method to retrieve previous agglomeration steps

  # Use feature names instead of indices
  if ( get.names ) {
    grouping <- lapply(grouping, function(x) {model@nodes[unlist(x)]})
  }

  # name the subnetworks
  names(grouping) <- paste("Subnet", 1:length(grouping), sep = "-")
  
  # If filters are given, apply them (stat needs to be specified)
  
  # SUBNET SIZE
  # treat size separately as it is quicker than responses
  
  if (any(c(!is.null(min.size), !is.null(max.size)))) {
    if ( is.null(min.size) ) { min.size <- 0   }
    if ( is.null(max.size) ) { max.size <- Inf }

    subnet.size <- sapply(grouping, length)    
    df <- data.frame(sapply(grouping, length))
    colnames(df) <- c("subnet.size")
    inds <- rownames(subset(df,
                            subnet.size >= min.size &
                            subnet.size <= max.size))

    grouping <- grouping[inds]

  }

  ## NUMBER OF RESPONSES
  if ( any(c( !is.null(min.responses)) ) ) {
    if ( is.null(stat) ) {
      cat("Computing subnetwork statistics, this may take a while..\n")
      stat <- result.stats(model) 
    }
    
    if ( is.null(min.responses) ) { min.responses <- 0 }

    # check which subnets pass the filter
    if ( length(inds) > 0 ) {
      inds <- rownames(subset(stat,
                              subnet.size >= min.size &
                              subnet.size <= max.size))#
                              #& subnet.responses >= min.responses))

      stat <- stat[inds, ]

      # Parse by number of responses
      stat <- stat[stat[, "subnet.responses"] >= min.responses,]
      
      
      grouping <- grouping[inds]
    } else { grouping <- NULL }
    
  }
  
  grouping
  
}

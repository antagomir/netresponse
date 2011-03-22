
setMethod("get.P.r", "NetResponseModel", function (model, subnet.id, log = TRUE) {

  # Prior probability for each response (mixture weights)
  # Pr(model, subnet.id)
  # output: a vector

  pars <- get.model.parameters(model, subnet.id)

  if (log) {
     log(pars$w)
   } else {
     pars$w
   }
})



setMethod("get.P.Sr", "NetResponseModel", function (sample, model, subnet.id, log = TRUE) {

  pars <- get.model.parameters(model, subnet.id)

  # Log density of a given sample group 
  # P(S|r) for each response; length of output equals to number of responses
  # i.e. logsum of the individual sample densities
  
  nodes <- model@subnets[[subnet.id]]

  # pick sample data for the response and 
  # ensure this is a matrix also when a single sample is given
  dat <- t(matrix(model@datamatrix[sample, nodes], ncol = length(nodes)))
  colnames(dat) <- sample
  rownames(dat) <- nodes
  
  P.Sr(dat, pars, log = log)
  
})


setMethod("get.P.rs.joint", "NetResponseModel", function (sample, model, subnet.id, log = TRUE) {
  
  # Joint probability P(r,s) where s can be a single point or set of
  # samples: P(r, s) = P(s | r) * P( r )
  pars <- get.model.parameters(model, subnet.id)  
  nodes <- model@subnets[[subnet.id]]
  # pick sample data for the response and 
  # ensure this is a matrix also when a single sample is given
  dat <- t(matrix(model@datamatrix[sample, nodes], ncol = length(nodes)))
  colnames(dat) <- sample
  rownames(dat) <- nodes

  P.rs.joint(dat, pars, log = TRUE)

})




setMethod("get.P.rS", "NetResponseModel", function (model, subnet.id, log = TRUE) {

  # Probability of a response, given sample (group)
  # P(r|S) = P(S|r)P(r)/P(S) = P(S, r)/(sum_r P(S, r))

  # Joint probability P(r,s) where s can be a single point or set of
  # samples: P(r, s) = P(s | r) * P( r )
  pars <- get.model.parameters(model, subnet.id)
  dat <- get.dat(model, subnet.id)
  
  P.rS(dat, pars, log = TRUE)

})


setMethod("get.P.rs.joint.individual", "NetResponseModel", function (sample, model, subnet.id, log = TRUE) {

  # Joint probability P(r,s) for individual samples
  # samples: P(r, s) = P(s | r) * P( r )
  pars <- get.model.parameters(model, subnet.id)  
  dat <- get.dat(model, subnet.id, sample)  # dat is now features x samples matrix
  
  prs <- P.rs.joint.individual(dat, pars, log = log)
  colnames(prs) <- sample

  prs
  
})




setMethod("get.P.s.individual", "NetResponseModel", function (sample, model, subnet.id, log = TRUE) {

  # Overall probability of sample s, given the model. 
  # individually for each sample

  # responses x samples
  # for each sample (column), density mass is the sum over joint densities on individual responses
  # P(s) = sum_r P(s, r) = sum_r P(s,r) = sum_r P(s|r)P(r)
  pars <- get.model.parameters(model, subnet.id)
  dat <- get.dat(model, subnet.id, sample)  # dat is now features x samples matrix
  
  ps <- P.s.individual(dat, pars, log = log)
  names(ps) <- sample

  ps
})




setMethod("sample.densities", "NetResponseModel", function (sample, model, subnet.id, log = TRUE, summarize = FALSE) {
    
  # Calculate conditional density  P(s|r) for each sample
  # in each response r and then in the complete model
  # if summarize = TRUE then give overall density of the sample,
  # otherwise densities for individual samples
  pars <- get.model.parameters(model, subnet.id)
  dat <- get.dat(model, subnet.id, sample)  # dat is now features x samples matrix
  
  psr <- P.sr(dat, pars, log = TRUE)

  # Density for the sample given the overall model
  ps <- P.s.individual(dat, pars, log = TRUE)
  names(ps) <- sample
  #dtot <- get.P.s.individual(sample, model, subnet.id, log = TRUE)
  # this should be same as colSums(psr * pars$w)

  rnam <- rownames(psr)
  psr <- rbind(psr, ps)
  rownames(psr) <- c(rnam, "total.density")
  
  # psr is a responses x samples matrix
  if (summarize) { psr <- rowSums(psr) } 
  
  if (!log) { psr <- exp(psr) }

  psr

})



setMethod("get.P.s", "NetResponseModel", function (sample, model, subnet.id, log = TRUE) {
  
  # Overall probability of sample s, given the Gaussian mixture model
  # P(s) = sum_r P(s, r)
  # ps <- log(sum(get.P.rs.joint(sample, model, subnet.id, log = FALSE)))

  # FIXME: numerically does not hold tightly that P(S) = sumr P(S, r) = sumr P(S|r)P(r)
  # the latter equality holds, P(S) is problematic. Differences are not big for examples
  # I checked, but they are still notable. Check in more detail this one.
  #sum(get.P.rs.joint(s, model, pars = NULL, subnet.id, log = FALSE))
  #sum(get.P.Sr(s, model, pars = NULL, subnet.id, log = FALSE) * pars$w)

  # P(s) separately for each individual sample
  
  # product over individual sample densities (i.e. log sum)
  #psi <- get.P.s.individual(sample, model, subnet.id, log = TRUE)
  pars <- get.model.parameters(model, subnet.id)
  dat <- get.dat(model, subnet.id, sample)  # dat is now features x samples matrix
    
  P.S(dat, pars, log = log)
  
})
  


setMethod("get.P.rs", "NetResponseModel", function (model, subnet.id, log = FALSE) {

  # Probability of each response, given sample
  # samples x responses matrix, each row sums to unity
  prs <- sample2response(model, subnet.id)
  
  if (log) {
     log(prs)
  } else {
    prs
  }
})



  
################################################################




setMethod("get.qofz", "NetResponseModel", function (model, subnet.id, log = FALSE) {

  # Retrieve P(r|s) from the model, given data and model parameters
  pars <- get.model.parameters(model, subnet.id)
  dat  <- get.dat(model, subnet.id)  # dat is now features x samples matrix
  qofz <- P.r.s(dat, pars, log = log)
  rownames(qofz) <- model@samples
  colnames(qofz) <- paste("Response", 1:ncol(qofz), sep = "-")   

  qofz
  
})


setMethod("get.dat", "NetResponseModel", function (model, subnet.id, sample = NULL) {

  if (is.null(sample)) {
    if (!is.null(rownames(model@datamatrix))) {
      sample <- rownames(model@datamatrix)
    } else {
      sample <- 1:nrow(model@datamatrix)
    }  
  }

  nodes <- model@subnets[[subnet.id]]
  dat <- t(matrix(model@datamatrix[sample, nodes], ncol = length(nodes)))
  rownames(dat) <- nodes
  colnames(dat) <- sample

  dat
})


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
 

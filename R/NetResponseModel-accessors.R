setMethod("get.P.r", "NetResponseModel", function(model, subnet.id, log = TRUE) {
    
    # Prior probability for each response (mixture weights) Pr(model, subnet.id)
    # output: a vector
    
    pars <- get.model.parameters(model, subnet.id)
    
    if (log) {
        log(pars$w)
    } else {
        pars$w
    }
})



setMethod("get.P.Sr", "NetResponseModel", function(sample, model, subnet.id, log = TRUE) {
    
    pars <- get.model.parameters(model, subnet.id)
    
    # Log density of a given sample group P(S|r) for each response; length of output
    # equals to number of responses i.e. logsum of the individual sample densities
    
    nodes <- model@subnets[[subnet.id]]
    
    # pick sample data for the response and ensure this is a matrix also when a
    # single sample is given
    dat <- t(matrix(model@datamatrix[sample, nodes], ncol = length(nodes)))
    colnames(dat) <- sample
    rownames(dat) <- nodes
    
    P.Sr(dat, pars, log = log)
    
})


setMethod("get.P.rs.joint", "NetResponseModel", function(sample, model, subnet.id, 
    log = TRUE) {
    
    # Joint probability P(r,s) where s can be a single point or set of samples: P(r,
    # s) = P(s | r) * P( r )
    pars <- get.model.parameters(model, subnet.id)
    nodes <- model@subnets[[subnet.id]]
    # pick sample data for the response and ensure this is a matrix also when a
    # single sample is given
    dat <- t(matrix(model@datamatrix[sample, nodes], ncol = length(nodes)))
    colnames(dat) <- sample
    rownames(dat) <- nodes
    
    P.rs.joint(dat, pars, log = TRUE)
    
})




setMethod("get.P.rS", "NetResponseModel", function(model, subnet.id, log = TRUE) {
    
    # Probability of a response, given sample (group) P(r|S) = P(S|r)P(r)/P(S) = P(S,
    # r)/(sum_r P(S, r))
    
    # Joint probability P(r,s) where s can be a single point or set of samples: P(r,
    # s) = P(s | r) * P( r )
    pars <- get.model.parameters(model, subnet.id)
    dat <- get.dat(model, subnet.id)
    
    P.rS(dat, pars, log = TRUE)
    
})


setMethod("get.P.rs.joint.individual", "NetResponseModel", function(sample, model, 
    subnet.id, log = TRUE) {
    
    # Joint probability P(r,s) for individual samples samples: P(r, s) = P(s | r) *
    # P( r )
    pars <- get.model.parameters(model, subnet.id)
    dat <- get.dat(model, subnet.id, sample)  # dat is now features x samples matrix
    
    prs <- P.rs.joint.individual(dat, pars, log = log)
    colnames(prs) <- sample
    
    prs
    
})




setMethod("get.P.s.individual", "NetResponseModel", function(sample, model, subnet.id, 
    log = TRUE) {
    
    # Overall probability of sample s, given the model.  individually for each sample
    
    # responses x samples for each sample (column), density mass is the sum over
    # joint densities on individual responses P(s) = sum_r P(s, r) = sum_r P(s,r) =
    # sum_r P(s|r)P(r)
    pars <- get.model.parameters(model, subnet.id)
    dat <- get.dat(model, subnet.id, sample)  # dat is now features x samples matrix
    
    ps <- P.s.individual(dat, pars, log = log)
    names(ps) <- sample
    
    ps
})




setMethod("sample.densities", "NetResponseModel", function(sample, model, subnet.id, 
    log = TRUE, summarize = FALSE) {
    
    # Calculate conditional density P(s|r) for each sample in each response r and
    # then in the complete model if summarize = TRUE then give overall density of the
    # sample, otherwise densities for individual samples
    pars <- get.model.parameters(model, subnet.id)
    dat <- get.dat(model, subnet.id, sample)  # dat is now features x samples matrix
    
    psr <- P.s.r(dat, pars, log = TRUE)
    
    # Density for the sample given the overall model
    ps <- P.s.individual(dat, pars, log = TRUE)
    names(ps) <- sample
    # dtot <- get.P.s.individual(sample, model, subnet.id, log = TRUE) this should be
    # same as colSums(psr * pars$w)
    
    rnam <- rownames(psr)
    psr <- rbind(psr, ps)
    rownames(psr) <- c(rnam, "total.density")
    
    # psr is a responses x samples matrix
    if (summarize) {
        psr <- rowSums(psr)
    }
    
    if (!log) {
        psr <- exp(psr)
    }
    
    psr
    
})



setMethod("get.P.s", "NetResponseModel", function(sample, model, subnet.id, log = TRUE) {
    
    # Overall probability of sample s, given the Gaussian mixture model P(s) = sum_r
    # P(s, r) ps <- log(sum(get.P.rs.joint(sample, model, subnet.id, log = FALSE)))
    
    # FIXME: numerically does not hold tightly that P(S) = sumr P(S, r) = sumr
    # P(S|r)P(r) the latter equality holds, P(S) is problematic. Differences are not
    # big for examples I checked, but they are still notable. Check in more detail
    # this one. sum(get.P.rs.joint(s, model, pars = NULL, subnet.id, log = FALSE))
    # sum(get.P.Sr(s, model, pars = NULL, subnet.id, log = FALSE) * pars$w)
    
    # P(s) separately for each individual sample
    
    # product over individual sample densities (i.e. log sum) psi <-
    # get.P.s.individual(sample, model, subnet.id, log = TRUE)
    pars <- get.model.parameters(model, subnet.id)
    dat <- get.dat(model, subnet.id, sample)  # dat is now features x samples matrix
    
    P.S(dat, pars, log = log)
    
})



setMethod("get.P.rs", "NetResponseModel", function(model, subnet.id, log = FALSE) {
    
    # Probability of each response, given sample samples x responses matrix, each row
    # sums to unity
    prs <- sample2response(model, subnet.id)
    
    if (log) {
        log(prs)
    } else {
        prs
    }
})


#' Sample-to-response matrix of probabilities P(r|s).
#' 
#' Retrieve P(r|s) from NetResponseModel model.
#' 
#' Calculates probability density for each response on a given sample based on
#' the estimated Gaussian mixture model.
#' 
#' @aliases getqofz getqofz,NetResponseModel-method
#' @usage getqofz(model, subnet.id, log = FALSE)
#' @param model NetResponseModel object.
#' @param subnet.id Subnetwork to investigate.
#' @param log Output in log probabilities.
#' @return Samples x responses matrix. Each entry is a probability P(r|s).
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation('netresponse').
#' @keywords utilities internal
#' @examples
#' 
#' # qofz <- getqofz(model, subnet.id, log = FALSE)
#' 
setMethod("getqofz", "NetResponseModel", function(model, subnet.id, log = FALSE) {
    
    # Retrieve P(r|s) from the model, given data and model parameters
    pars <- get.model.parameters(model, subnet.id)
    dat <- get.dat(model, subnet.id)  # Dat is now features x samples matrix
    qofz <- P.r.s(dat, pars, log = log)
    rownames(qofz) <- rownames(model@datamatrix)  #model@samples
    colnames(qofz) <- paste("Mode", seq_len(ncol(qofz)), sep = "-")
    
    qofz
    
})


#' Get subnetwork data
#' 
#' @inheritParams sample2response
#' @param sample Define the retrieved samples
#' @aliases get.dat get.dat,NetResponseModel-method
#' @return Subnet data matrix
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation('netresponse')
#' @export
#' @keywords utilities
#' @examples
#' ## Load a pre-calculated netresponse model obtained with 
#' # model <- detect.responses(toydata$emat, toydata$netw, verbose = FALSE)
#' # data( toydata ); get.dat(toydata$model) 
setMethod("get.dat", "NetResponseModel", function(model, subnet.id, sample = NULL) {
    
    # usage get.dat(model, subnet.id, sample = NULL)
    
    if (is.null(sample)) {
        sample <- seq_len(nrow(model@datamatrix))
    }
    nodes <- model@subnets[[subnet.id]]
    
    dat <- t(matrix(model@datamatrix[sample, nodes], length(sample)))
    rownames(dat) <- nodes
    colnames(dat) <- as.character(sample)
    
    dat
})


#' get.subnets
#' 
#' List the detected subnetworks (each is a list of nodes in the corresponding
#' subnetwork).
#' 
#' @aliases get.subnets get.subnets,NetResponseModel-method
#' @param model Output from the detect.responses function. An object of
#' NetResponseModel class.
#' @param get.names Logical. Indicate whether to return subnetwork nodes using
#' node names (TRUE) or node indices (FALSE).
#' @param min.size,max.size Numeric. Filter out subnetworks whose size is not
#' within the limits specified here.
#' @param min.responses Numeric. Filter out subnetworks with less responses
#' (mixture components) than specified here.
#' @return A list of subnetworks.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references Leo Lahti et al.: Global modeling of transcriptional responses
#' in interaction networks. Bioinformatics (2010).  See citation('netresponse')
#' for details.
#' @keywords utilities
#' @export
#' @examples
#' ## Load a pre-calculated netresponse model obtained with 
#' # model <- detect.responses(toydata$emat, toydata$netw, verbose = FALSE)
#' # data( toydata ); get.subnets(toydata$model) 
setMethod("get.subnets", "NetResponseModel", function(model, get.names = TRUE, min.size = 2, 
    max.size = Inf, min.responses = 2) {
    
    grouping <- model@last.grouping
    
    # Use feature names instead of indices
    if (get.names) {
        grouping <- lapply(grouping, function(x) {
            colnames(model@datamatrix)[unlist(x)]
        })
    }
    
    # If filters are given, apply them (stat needs to be specified)
    
    # SUBNET SIZE
    
    subnet.size <- sapply(grouping, length)
    df <- data.frame(sapply(grouping, length))
    colnames(df) <- c("subnet.size")
    inds <- rownames(subset(df, subnet.size >= min.size & subnet.size <= max.size))
    
    ## NUMBER OF RESPONSES
    
    if (min.responses > 1) {
        
        stat <- model.stats(model)
        
        # check which subnets pass the filter
        inds.size <- rownames(subset(stat, subnet.size >= min.size & subnet.size <= 
            max.size))  # & 
        
        inds.nresp <- rownames(stat)[which(stat[["responses"]] >= min.responses)]
        
        inds <- intersect(inds.size, inds.nresp)
        
    }
    
    # Get the filtered subnetwork list
    if (length(inds) == 0) {
        grouping <- NULL
    } else {
        grouping <- grouping[inds]
    }
    
    grouping
    
})


setMethod(f = "[[", signature("NetResponseModel"), definition = (function(x, i, j = "missing", 
    ..., exact = TRUE) {
    if (typeof(i) == "numeric") {
        i <- names(x)[[i]]
    }
    get.model.parameters(x, subnet.id = i)
}))

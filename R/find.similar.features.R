find.similar.features <- function (model, subnet.id, datamatrix = NULL, verbose = FALSE) {

  # Given subnetwork, order the remaining genes in the data by
  # similarity with this subnetwork (coordinated transcritional response)
  # the same similarity score is used here than in subnetwork agglomeration steps

  # By default, investigate the same data matrix which was used for modelling
  # NOTE: default parameters of the model are always used (model@params)
  if (is.null(datamatrix)) {
     datamatrix <- model@datamatrix
  }

  subnetwork.nodes <- model@subnets[[subnet.id]]
  other.nodes <- setdiff(colnames(datamatrix), subnetwork.nodes)

  # Check that subnetwork nodes are included in the data
  if (!all(subnetwork.nodes %in% colnames(datamatrix))) {
    stop("The complete subnetwork (features) must be included in the input datamatrix (samples x features)")
  }

  # Use indices of the nodes instead of names
  other.nodes.idx <- match(other.nodes, colnames(datamatrix))
  subnet.idx <- match(subnetwork.nodes, colnames(datamatrix))
  Nlog <- log(nrow(datamatrix)) # log of sample size
    
  #############################################################

  # FIXME: take straight from the model
  # Retrieve model for the investigated subnetwork
  m.subnet <- vdp.mixt( matrix(datamatrix[, subnetwork.nodes], nrow( datamatrix )),
                  implicit.noise = model@params$implicit.noise,
		  prior.alpha = model@params$prior.alpha,
                  prior.alphaKsi = model@params$prior.alphaKsi,
		  prior.betaKsi = model@params$prior.betaKsi,
		  threshold = model@params$threshold, initial.K =
		  model@params$initial.responses,
                  ite = model@params$ite, 
		  c.max = model@params$max.responses - 1 )
                  
  bic.subnet  <- m.subnet$posterior$Nparams*Nlog + 2*m.subnet$free.energy # BIC for this model
     
  ############################################################

  # Go through data points (features i.e. nodes) and measure similarity for each
  # feature with the given subnetwork
  delta <- c()
  for (fi in other.nodes.idx) {
    fname <- colnames(datamatrix)[[fi]] # pick also node name
    
    # Note: this calculation jumps the indices that are in the subnetwork
    if (verbose) { message(paste("Calculating similarity for feature ", fi, "/", max(other.nodes.idx))) }
  
    #print(" Retrieve models for the individual node fi")
        
     m.node <- vdp.mixt( matrix(datamatrix[, fi], nrow( datamatrix )),
                  implicit.noise = model@params$implicit.noise, prior.alpha = model@params$prior.alpha,
                  prior.alphaKsi = model@params$prior.alphaKsi, prior.betaKsi = model@params$prior.betaKsi,
                  threshold = model@params$threshold, initial.K = model@params$initial.responses,
                  ite = model@params$ite, c.max = model@params$max.responses - 1 )

     bic.node <- m.node$posterior$Nparams*Nlog + 2*m.node$free.energy 
     # FIXME: make a function to readily give BIC (or other objective function for a particular model?)

    ##########################################################################################

    #print(" JOINT MODEL")

    # subnetwork + candidate feature
    vars  <- sort(c(subnet.idx, fi))

    m.joint <- vdp.mixt(matrix(datamatrix[, vars], nrow( datamatrix )),
                    implicit.noise = model@params$implicit.noise, prior.alpha = model@params$prior.alpha,
                    prior.alphaKsi = model@params$prior.alphaKsi, prior.betaKsi = model@params$prior.betaKsi,
                    threshold = model@params$threshold, initial.K = model@params$initial.responses,
                    ite = model@params$ite, c.max = model@params$max.responses - 1 )

    bic.joint  <- m.joint$posterior$Nparams*Nlog + 2*m.joint$free.energy
    # = combined: subnet + the new node

    #####################################################################################################

    #print(" TWO INDEPENDENT MODELS")
    bic.ind   <- bic.subnet + bic.node

    #####################################################################################

    #print(" BIC-cost for two independent vs. joint model")
    # change (increase) of the total BIC (cost)
    delta[[fname]] <- bic.joint - bic.ind

  }

  # Return the BIC delta values. The smaller this is, the more similar 
  # the feature is in regard to the subnetwork
  # negative values mean that the subnetwork model would be improved by merging the gene (= feature)
  df <- as.data.frame(delta)
  df[["feature.name"]] <- rownames(df)
  df <- df[c("feature.name", "delta")]
  df <- fsort(df, "delta")
  
  df
}




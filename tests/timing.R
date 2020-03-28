
# Play with different options and check their effect on  running times for bic and vdp 

skip <- TRUE

if (!skip) {

  Ns <- 100
  Nd <- 2

  set.seed(3488400)

  D <- cbind(

     	rbind(matrix(rnorm(Ns*Nd, mean = 0), ncol = Nd), 
       	      matrix(rnorm(Ns*Nd, mean = 2), ncol = Nd),
      	      cbind(rnorm(Ns, mean = -1), rnorm(Ns, mean = 3))
 	    ), 

     	rbind(matrix(rnorm(Ns*Nd, mean = 0), ncol = Nd), 
       	      matrix(rnorm(Ns*Nd, mean = 2), ncol = Nd),
      	      cbind(rnorm(Ns, mean = -1), rnorm(Ns, mean = 3))
 	    )
	    )

  rownames(D) <- paste("R", 1:nrow(D), sep = "-")
  colnames(D) <- paste("C", 1:ncol(D), sep = "-")

  ts <- c()
  for (mm in c("bic", "vdp")) {


    # NOTE: no PCA basis needed with mixture.method = "bic"
    tt <- system.time(detect.responses(D, verbose = TRUE, max.responses = 5, 
	   		       mixture.method = mm, information.criterion = "BIC", 
			       merging.threshold = 0, bic.threshold = 0, pca.basis = TRUE))

    print(paste(mm, ":", round(tt[["elapsed"]], 3)))
    ts[[mm]] <- tt[["elapsed"]]
  }

   print(paste(names(ts)[[1]], "/", names(ts)[[2]], ": ", round(ts[[1]]/ts[[2]], 3)))

}

# -> VDP is much faster when sample sizes increase 
# 1000 samples -> 25-fold speedup with VDP



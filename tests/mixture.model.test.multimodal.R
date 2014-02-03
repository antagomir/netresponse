library(netresponse)

# Three MODES

# set.seed(34884)
set.seed(3488400)

Ns <- 200
Nd <- 2

D3 <- rbind(matrix(rnorm(Ns*Nd, mean = 0), ncol = Nd), 
      	    matrix(rnorm(Ns*Nd, mean = 3), ncol = Nd),
      	    cbind(rnorm(Ns, mean = -3), rnorm(Ns, mean = 3))
	    )

X11()
par(mfrow = c(2,2))
for (mm in c("vdp", "bic")) {
  for (pp in c(FALSE, TRUE)) {

    # Fit nonparametric Gaussian mixture model
    out <- mixture.model(D3, mixture.method = mm, pca.basis = pp)
    plot(D3, col = apply(out$qofz, 1, which.max), main = paste(mm, "/ pca:",  pp)) 

  }
}

# VDP is less sensitive than BIC in detecting Gaussian modes (more
# separation between the clusters needed)

# pca.basis option is less important for sensitive detection but
# it will help to avoid overfitting to unimodal features that
# are not parallel to the axes (unimodal distribution often becomes
# splitted in two or more clusters in these cases)


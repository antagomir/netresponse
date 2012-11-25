library(netresponse)

# SINGLE MODE

# Produce test data that has full covariance
# It is expected that
# pca.basis = FALSE splits Gaussian with full covariance into two modes
# pca.basis = TRUE should detect just a single mode

Ns <- 200
Nd <- 2
k <- 1.5

D2 <- matrix(rnorm(Ns*Nd), ncol = Nd) %*% rbind(c(1,k), c(k,1))

par(mfrow = c(2,2))
for (mm in c("vdp", "bic")) {
  for (pp in c(FALSE, TRUE)) {

    # Fit nonparametric Gaussian mixture model
    out <- mixture.model(D2, mixture.method = mm, pca.basis = pp)
    plot(D2, col = apply(out$qofz, 1, which.max), main = paste("mm:" , mm, "/ pp:",  pp)) 

  }
}


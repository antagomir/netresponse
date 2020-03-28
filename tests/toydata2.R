# Generate Nc components from normal-inverseGamma prior

set.seed(12346)

Ns <- 300
Nd <- 2

# Isotropic cloud
D1 <- matrix(rnorm(Ns*Nd), ncol = Nd) 

# Single diagonal mode
D2 <- matrix(rnorm(Ns*Nd), ncol = Nd) %*% rbind(c(1,2), c(2,1)) 

# Two isotropic modes
D3 <- rbind(matrix(rnorm(Ns/2*Nd), ncol = Nd), matrix(rnorm(Ns/2*Nd, mean = 3), ncol = Nd))
D <- cbind(D1, D2, D3)

colnames(D) <- paste("Feature-",  1:ncol(D), sep = "")
rownames(D) <- paste("Sample-", 1:nrow(D), sep = "")


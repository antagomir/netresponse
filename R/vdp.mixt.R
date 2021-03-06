# system('~/local/R/R-2.13.0/bin/R CMD SHLIB ../src/netresponse.c');
# dyn.load('../src/netresponse.so')

#' vdp.mixt
#' 
#' Accelerated variational Dirichlet process Gaussian mixture.
#' 
#' Implementation of the Accelerated variational Dirichlet process Gaussian
#' mixture model algorithm by Kenichi Kurihara et al., 2007.
#'
#' ALGORITHM SUMMARY 
#' This code implements Gaussian mixture models with diagonal covariance matrices. 
#' The following greedy iterative approach is taken in order to obtain the number
#' of mixture models and their corresponding parameters:
#'
#' 1. Start from one cluster, $T = 1$.
#' 2. Select a number of candidate clusters according to their values of 
#'    'Nc' = \\sum_{n=1}^N q_{z_n} (z_n = c) (larger is better).
#' 3. For each of the candidate clusters, c: 
#'     3a. Split c into two clusters, c1 and c2, through the bisector of its 
#'         principal component. Initialise the responsibilities 
#'         q_{z_n}(z_n = c_1) and q_{z_n}(z_n = c_2). 
#'     3b. Update only the parameters of c1 and c2 using the observations that
#'         belonged to c, and determine the new value for the free energy, F{T+1}.
#'     3c. Reassign cluster labels so that cluster 1 corresponds to the largest 
#'         cluster, cluster 2 to the second largest, and so on.
#' 4. Select the split that lead to the maximal reduction of free energy, F{T+1}.
#' 5. Update the posterior using the newly split data.
#' 6. If FT - F{T+1} < \\epsilon then halt, else set T := T +1 and go to step 2.
#'
#' The loop is implemented in the function greedy(...)
#'
#' @param dat Data matrix (samples x features).
#' @param prior.alpha,prior.alphaKsi,prior.betaKsi Prior parameters for
#' Gaussian mixture model (normal-inverse-Gamma prior). alpha tunes the mean;
#' alphaKsi and betaKsi are the shape and scale parameters of the inverse Gamma
#' function, respectively.
#' @param do.sort When true, qOFz will be sorted in decreasing fashion by
#' component size, based on colSums(qOFz). The qOFz matrix describes the
#' sample-component assigments in the mixture model.
#' @param threshold Defines the minimal free energy improvement that stops the
#' algorithm: used to define convergence limit.
#' @param initial.K Initial number of mixture components.
#' @param ite Defines maximum number of iterations on posterior update
#' (updatePosterior). Increasing this can potentially lead to more accurate
#' results, but computation may take longer.
#' @param implicit.noise Adds implicit noise; used by vdp.mk.log.lambda.so and
#' vdp.mk.hp.posterior.so. By adding noise (positive values), one can avoid
#' overfitting to local optima in some cases, if this happens to be a problem.
#' @param c.max Maximum number of candidates to consider in
#' find.best.splitting. During mixture model calculations new mixture
#' components can be created until this upper limit has been reached. Defines
#' the level of truncation for a truncated stick-breaking process.
#' @param speedup When learning the number of components, each component is
#' splitted based on its first PCA component. To speed up, approximate by using
#' only subset of data to calculate PCA.
#' @param min.size Minimum size for a component required for potential
#'   splitting during mixture estimation.
#'
#' @return \item{ prior }{Prior parameters of the vdp-gm model (qofz: priors on observation lables; Mu: centroids; S2: variance).} 
#'         \item{ posterior }{Posterior estimates for the model parameters and statistics.} 
#'         \item{ weights }{Mixture proportions, or weights, for the Gaussian mixture components.}
#'  \item{ centroids }{Centroids of the mixture components.} 
#'  \item{ sds }{ Standard deviations for the mixture model components (posterior modes of the covariance diagonals square root). Calculated as sqrt(invgam.scale/(invgam.shape + 1)). } 
#'  \item{ qOFz }{ Sample-to-cluster assigments (soft probabilistic associations).} 
#' \item{ Nc }{Component sizes}
#'  \item{ invgam.shape }{ Shape parameter (alpha) of the inverse Gamma distribution } 
#'  \item{ invgam.scale }{ Scale parameter (beta) of the inverse Gamma distribution } 
#'  \item{ Nparams }{ Number of model parameters }
#' \item{ K }{ Number of components in the mixture model } 
#'  \item{ opts }{Model parameters that were used.} 
#'  \item{ free.energy }{Free energy of the model.}
#'
#' @note This implementation is based on the Variational Dirichlet Process
#'   Gaussian Mixture Model implementation, Copyright (C) 2007 Kenichi Kurihara
#'   (all rights reserved) and the Agglomerative Independent Variable Group
#'   Analysis package (in Matlab): Copyright (C) 2001-2007 Esa Alhoniemi, Antti
#'   Honkela, Krista Lagus, Jeremias Seppa, Harri Valpola, and Paul Wagner.
#' @author Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references Kenichi Kurihara, Max Welling and Nikos Vlassis: Accelerated
#'   Variational Dirichlet Process Mixtures.  In B. Sch\'olkopf and J. Platt and
#'   T. Hoffman (eds.), Advances in Neural Information Processing Systems 19,
#'   761--768. MIT Press, Cambridge, MA 2007.
#' @keywords methods iteration
#' @export
#' @examples
#' 
#'   set.seed(123)
#' 
#'   # Generate toy data with two Gaussian components
#'   dat <- rbind(array(rnorm(400), dim = c(200,2)) + 5,
#'                array(rnorm(400), dim = c(200,2)))
#' 
#'   # Infinite Gaussian mixture model with 
#'   # Variational Dirichlet Process approximation
#'   mixt <- vdp.mixt( dat )
#' 
#'   # Centroids of the detected Gaussian components
#'   mixt$posterior$centroids
#' 
#'   # Hard mixture component assignments for the samples
#'   apply(mixt$posterior$qOFz, 1, which.max)
#' 
#' 
vdp.mixt <- function(dat, prior.alpha = 1, prior.alphaKsi = 0.01, prior.betaKsi = 0.01, 
    do.sort = TRUE, threshold = 1e-05, initial.K = 1, ite = Inf, implicit.noise = 0, 
    c.max = 10, speedup = TRUE, min.size = 5) {
    
    # qOFz sorted in decreasing fashion based on colSums(qOFz) threshold: minimal
    # free energy improvement that stops the algorithm initial.K initial number of
    # components ite # used on updatePosterior: maximum number of iterations
    # implicit.noise # Adds implicit noise in vdp.mk.log.lambda.so and
    # vdp.mk.hp.posterior.so
    
    # c.max max. candidates to consider in Candidates are chosen based on their Nc
    # value (larger = better). Nc = colSums(qOFz) find.best.splitting. i.e.
    # truncation parameter speedup: during DP, components are splitted based on their
    # first PCA component. To speed up, approximate by using only subset data to
    # calculate PCA.  min.size # min size for a component to be splitted
    
    # Prior parameters
    opts <- list(prior.alpha = prior.alpha, prior.alphaKsi = prior.alphaKsi, prior.betaKsi = prior.betaKsi, 
        speedup = speedup, do.sort = do.sort, threshold = threshold, initial.K = initial.K, 
        ite = ite, implicitnoisevar = implicit.noise, c.max = c.max)
    
    # alphaKsi smaller -> less clusters (and big!) -> quite sensitive betaKsi larger
    # -> less clusters (and big!)
    
    data <- list()
    data[["given.data"]] <- list()
    data[["given.data"]][["X1"]] <- dat
    
    # sample-component assignments
    qOFz <- rand.qOFz(nrow(dat), initial.K)
    
    # The hyperparameters of priors
    hp.prior <- mk.hp.prior(data, opts)
    
    # Posterior
    hp.posterior <- mk.hp.posterior(data, qOFz, hp.prior, opts)
    
    # Note: greedy gives components in decreasing order by size
    templist <- greedy(data, hp.posterior, hp.prior, opts, min.size)
    templist$hp.prior <- c(templist$hp.prior, list(qOFz = qOFz))
    qOFz <- matrix(templist$hp.posterior$qOFz, nrow(dat))
    
    ############################################### 
    
    # Retrieve model parameters
    
    # number of mixture components (nonempty components only!)  response must have
    # non-negligible' probability mass!  i.e. at least some points associated with it
    Kreal <- max(apply(qOFz, 1, which.max))  #sum(colSums(qOFz) > 1e-3)
    qOFz <- matrix(qOFz[, seq_len(Kreal)], nrow(dat))
    rownames(qOFz) <- rownames(dat)
    
    # Calculate mixture model parameters FIXME: move this outside from this vdp.mixt
    # function
    
    # Negative free energy is (variational) lower bound for P(D|H) Use it to
    # approximate P(D|HClist <- list(C))
    
    # number of parameters (d-dim. centroid + diag. cov. matrix + component weight
    # for each Kreal)
    Nparams <- Kreal * (2 * ncol(dat) + 1)
    
    # Calculate map estimates of model parameters from the posterior
    
    # variances are assumed inverse Gamma distributed and here beta/alpha gives the
    # expectation)
    
    # Parameters of the inverse Gamma function for component variances
    invgam.shape <- matrix(templist$hp.posterior$KsiAlpha[seq_len(Kreal), ], Kreal)
    invgam.scale <- matrix(templist$hp.posterior$KsiBeta[seq_len(Kreal), ], Kreal)
    
    # Calculate variances (mean and mode of the invgam distr.) from scale and shape
    # FIXME: beta/alpha used in C code var.update <-
    # matrix(invgam.scale/invgam.shape, Kreal) var.mean <-
    # matrix(invgam.scale/(invgam.shape - 1), Kreal)
    var.mode <- matrix(invgam.scale/(invgam.shape + 1), Kreal)
    variances <- var.mode  # select mean, mode, or their average for update
    
    # Ignore empty components assuming that the components have been ordered in
    # decreasing order by size component centroids
    centroids <- matrix(templist$hp.posterior$Mubar[seq_len(Kreal), ], Kreal)
    
    ############################################# 
    
    # Estimate weights directly from mass on each component; with large sample size
    # the solution should converge there
    w <- colSums(qOFz)
    w <- w/sum(w)  # normalize
    # FIXME: improve later, or retrieve weight directly from vdp code or the
    # densities to make more robust also for small sample size NOTE: the way we
    # calculate weights here is not in exact agreement with the real weights in the
    # model that would give qOFz but it is an approximation, and needed for end
    # analysis
    
    posterior <- list(weights = w, centroids = centroids, sds = sqrt(variances), 
        qOFz = qOFz, Nc = colSums(qOFz), invgam.shape = invgam.shape, invgam.scale = invgam.scale, 
        Nparams = Nparams, K = Kreal)
    
    # centroids: mubar sds alpha, beta Nc: # component sizes invgam shape # KsiAlpha
    # invgam.scale # KsiBeta Nparams = Nparams, # number of model parameters K =
    # Kreal # number of components
    
    # Later: include these from hp.posterior to the output later if needed.
    # 'Mutilde[1:Kreal,]' 'gamma[,1:Kreal]' 'Uhat'
    
    list(prior = templist$hp.prior, posterior = posterior, opts = opts, free.energy = templist$free.energy)
}

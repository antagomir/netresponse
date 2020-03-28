#' @title pick.model.parameters
#' @description Pick model parameters
#' @param m vdp.mixt output
#' @param nodes node names for naming purposes 
#' @return Model parameters
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse")
#' @keywords internal
#' @export
#' @examples #
pick.model.parameters <- function (m, nodes) {
 
  # Pick parameters
  w    <- m$posterior$weights    # component weights
  mu   <- m$posterior$centroids  # component centroids
  sds  <- m$posterior$sds        # component standard devs
  qofz  <- m$posterior$qOFz      # soft mode assignmentd
  free.energy <- m$free.energy # free energy
  Nparams <- m$posterior$Nparams

  # For mu and std, rows correspond to the mixture components, in w the elements
  list(mu = mu, sd = sds, w = w, free.energy = free.energy, Nparams = Nparams, qofz = qofz)

}



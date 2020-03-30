
#' Convert grouping info into a list; each element corresponds to a
#' group and lists samples in that group. 
#' 
#' @param groupings a list, a vector, or a samplesxmodes assignment matrix
#' @param verbose verbose
#'
#' @return Group list
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation('netresponse')
#' @keywords utilities
#' @export
#' @examples res <- listify.groupings(groupings = NULL)
listify.groupings <- function(groupings, verbose = FALSE) {

    if (is.null(groupings)) {return(NULL)}

    if (is.matrix(groupings[[1]])) {
        if (verbose) {
            message("Convert mode matrix to list")
        }
        groupings.list <- list()
        for (sn in names(groupings)) {
            # samples in each mode (hard assignment)
            groupings[[sn]] <- response2sample(groupings[[sn]])
        }
    } else if (is.vector(groupings) && !is.list(groupings)) {
        if (verbose) {
            message("Convert mode vector to list")
        }
        groupings.list <- split(names(groupings), groupings)
    } else if (is.list(groupings)) {
        if (verbose) {
            message("Mode list ok")
        }
        groupings.list <- groupings
    }
    
    groupings.list
}


#' Convert grouping info into a vector; each element corresponds to a
#' group and lists samples in that group.
#' 
#' @param groupings a list, a vector, or a samplesxmodes assignment matrix
#' @param verbose verbose
#'
#' @return Indicator vector
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation('netresponse')
#' @keywords utilities
#' @export
#' @examples #

vectorize.groupings <- function(groupings, verbose = FALSE) {
    
    if (is.matrix(groupings[[1]])) {
        if (verbose) {
            message("Convert mode matrix to vector")
        }
        gv <- colnames(groupings)[apply(groupings, 1, which.max)]
    } else if (is.vector(groupings) && !is.list(groupings)) {
        if (verbose) {
            message("Convert mode vector to list")
        }
        gv <- groupings
    } else if (is.list(groupings)) {
        gv <- unlist(sapply(seq_len(length(groupings)), function(i) {
            rep(names(groupings)[[i]], length(groupings[[i]]))
        }))
        names(gv) <- unlist(groupings)
    }
    
    gv
}

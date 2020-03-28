#' @title enrichment.list.factor
#' @description Orders the responses by association strength (enrichment score) to a given sample set. For instance, if the samples correspond to a particular experimental factor, this function can be used to prioritize the responses according to their association strength to this factor.
#' @param groupings List of groupings. Each model should have a sample-cluster assignment matrix qofz.
#' @param method 'hypergeometric' measures enrichment of factor levels in this
#'   response; 'precision' measures response purity for each factor level;
#'   'dependency' measures logarithm of the joint density between response and
#'   factor level vs. their marginal densities: log(P(r,s)/(P(r)P(s)))
#' @param verbose Follow progress by intermediate messages.
#' @param annotation.vector annotation vector
#' @param level level
#' @return A data frame which gives a data
#'   frame of responses ordered by enrichment score for the investigated sample.
#'   The model, response id and enrichment score are shown. The method field
#'   indicates the enrichment calculation method. The sample field lists the
#'   samples et for which the enrichments were calculated. The info field lists
#'   additional information on enrichment statistics.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse") for citation details.
#' @keywords utilities
#' @importFrom qvalue qvalue
#' @export
#' @examples #
enrichment.list.factor.minimal <- function (groupings, method, verbose = FALSE, annotation.vector, level) {

  # groupings, level.samples, method = method
  level.samples <- names(annotation.vector)[annotation.vector %in% level]

  if (is.null(names(groupings))) { names(groupings) <- 1:length(groupings) }

  # Check enrichment in the selected responses  
  enrichment.info <- list()
  cnt <- 0
  for (subnet.id in names(groupings)) {

    groupings.list <- listify.groupings(groupings[[subnet.id]])

    for (response in names(groupings.list)) {

      # In enrichment calculations use only those samples where annotations are available
      total.samples <- intersect(unlist(groupings.list), names(which(!is.na(annotation.vector))))
      response.samples <- intersect(groupings.list[[response]], names(which(!is.na(annotation.vector))))

      enr <- response.enrichment(total.samples, response.samples, level.samples, method)
      
      cnt <- cnt + 1
      enrichment.info[[cnt]] <- c(model = subnet.id, mode = response, enrichment.score = enr$score, enr$info) 

    }
  }
    
  if (verbose) { message("Groupings checked.") }

  if (length(enrichment.info) > 0) {

    enrichment.info <- enrichment.info[sapply(enrichment.info, function (ei) {length(ei) > 2})]

    enr <- as.data.frame(t(sapply(enrichment.info, identity)))

    if ("pvalue" %in% colnames(enr)) {
        
      if (length(enr$pvalue) > 100) {
        # calculate q-values
        enr$qvalue <- qvalue(as.numeric(as.character(enr$pvalue)))$qvalues
      } else if (length(enr$pvalue) > 10) {
        enr$qvalue <- qvalue(as.numeric(as.character(enr$pvalue)), pi0.method = "bootstrap", fdr.level = 0.25)$qvalues
      } else {

        warning("Not enough p-values for q-value estimation")
        enr$qvalue <- rep(NA, length(enr$pvalue))

      }
    }

    if (nrow(enr) == 0) { return(NULL) }

    enr[,3:ncol(enr)] <- apply(enr[,3:ncol(enr)], 2, as.numeric)

    if ("enrichment.score" %in% names(enr)) {

      enr <- enr[order(enr$enrichment.score, decreasing = TRUE),]

      enr[["model"]] <- as.character(enr[["model"]])

      tmp <- list(ordered.responses = enr, method = method, sample = level.samples)

      return(tmp)

    } else {
      return(NULL)
    }

  } else {
    return(NULL)
  }	

}


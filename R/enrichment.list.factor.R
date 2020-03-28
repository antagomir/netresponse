#' @title enrichment.list.factor
#' @description Orders the responses by association strength (enrichment score) to a given sample set. For instance, if the samples correspond to a particular experimental factor, this function can be used to prioritize the responses according to their association strength to this factor.
#' @param models List of models. Each model should have a sample-cluster assignment matrix qofz.
#' @param level.samples Measure enrichment of this sample (set) across the observed responses.
#' @param method 'hypergeometric' measures enrichment of factor levels in this
#'   response; 'precision' measures response purity for each factor level;
#'   'dependency' measures logarithm of the joint density between response and
#'   factor level vs. their marginal densities: log(P(r,s)/(P(r)P(s)))
#' @param verbose Follow progress by intermediate messages.
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
enrichment.list.factor <- function (models, level.samples, method, verbose = FALSE) {

  # models; level.samples <- level.samples; method = method		       		       
  if (is.null(names(models))) { names(models) <- 1:length(models) }

  # Get model statistics
  stat <- model.stats(models)
	       
  # Check enrichment in the selected responses  
  enrichment.info <- list()
  for (subnet.id in names(models)) {

    if ( verbose ) { message(subnet.id) }
    
    qofz <- models[[subnet.id]]$qofz

    for (response in 1:ncol(qofz)) {

      enr <- response.enrichment(qofz, level.samples, response, method)

      # add further info about enrichments
      cnt <- cnt + 1

      enrichment.info[[cnt]] <- c(model = subnet.id, mode = response, enrichment.score = enr$score, enr$info) 

    }
  }
    
  if (verbose) { message("Models checked.") }

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

    enr[,3:ncol(enr)] <- apply(enr[,3:ncol(enr)], 2, as.numeric)

    if ("enrichment.score" %in% names(enr)) {

      enr <- enr[order(enr$enrichment.score, decreasing = TRUE),]

      enr[["model"]] <- as.character(enr[["model"]])

      # Add subnet info in the result table
      enr <- cbind(enr, stat[enr$model,])

      tmp <- list(ordered.responses = enr, method = method, sample = level.samples)

      return(tmp)

    } else {
      return(NULL)
    }

  } else {
    return(NULL)
  }	

}


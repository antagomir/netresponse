# Fixme: finish this later
#order.samples <- function (subnet.id, model, phenodata, which.factor, response, method = "hypergeometric") {
#    
#  # - for given response, order factor levels by association strength (enrichment score)
#  #   P(s|r) = P(s,r)/P(r) total sample density and/or average sample density (both for individuals and groups s/S)
#  #   P(s|r)/P(s) = P(s,r)/P(r)P(s) 
#
#  #   P(s|r)/P(s) -> OK: method = "dependency"
#  #   P(s|r) = P(s,r)/P(r)  -> OK, normalized version P(s|r)/P(S|r) available: method = "precision"
#
#  enr <- response.enrichments(subnet.id, model, phenodata, which.factor, response, method)
# 
#  # For the given response, return levels of the given factor (decreasing ordering by enrichement score)
#  sort(enr, decreasing = TRUE)
#}


#' @title order.responses
#' @description Orders the responses by association strength (enrichment score) to a given sample set. For instance, if the samples correspond to a particular experimental factor, this function can be used to prioritize the responses according to their association strength to this factor.
#' @param models List of models. Each model should have a sample-cluster assignment matrix qofz.
#' @param sample Measure enrichment of this sample (set) across the observed
#'   responses.
#' @param method 'hypergeometric' measures enrichment of factor levels in this
#'   response; 'precision' measures response purity for each factor level;
#'   'dependency' measures logarithm of the joint density between response and
#'   factor level vs. their marginal densities: log(P(r,s)/(P(r)P(s)))
#' @param min.size,max.size,min.responses Optional parameters to filter the
#'   results based on subnet size and number of responses.
#' @param subnet.ids Specify subnets for which the responses shall be ordered.
#'   By default, use all subnets.
#' @param verbose Follow progress by intermediate messages.
#' @param data data (samples x features; or a vector in univariate case)
#' @return A data frame with elements 'ordered.responses' which gives a data
#'   frame of responses ordered by enrichment score for the investigated sample.
#'   The subnetwork, response id and enrichment score are shown. 
#'   The method field indicates the enrichment calculation method. 
#'   The sample field lists the samples et for which the enrichments 
#'   were calculated. 
#'   The info field lists additional information on enrichment statistics.
#' @note Tools for analyzing end results of the model.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @references See citation("netresponse") for citation details.
#' @keywords utilities
#' @importFrom qvalue qvalue
#' @export
#' @examples #
#' # - for given sample/s (factor level), 
#' #     order responses (across all subnets) by association strength 
#' #     (enrichment score); overrepresentation
#' # order.responses(model, sample, method  = "hypergeometric") 
order.responses <- function (models, sample, method = "hypergeometric", 
		   	     min.size = 2, max.size = Inf, 
			     min.responses = 2, 
			     subnet.ids = NULL, 
			     verbose = FALSE, data = NULL) {

    # Given sample (for instance set of samples associated with a given factor 
    # level) order the responses across all subnetworks based on their 
    # association strength

    # For a given sample, calculate enrichment values in each response
    subnets <- responses <- scores <- pvals <- c()
    enrichment.info <- list()
    cnt <- 0

    # Get model statistics
    stat <- model.stats(models)

    # Filter the results
    sn <- get.subnets(models, get.names = TRUE, min.size, max.size, min.responses)
    stat <- stat[names(sn),]
    if (is.null(subnet.ids)) { subnet.ids <- rownames(stat) }

    enr <- enrichment.list.factor(models, sample, method)

  enr

}
			    
#' @title Listing significant responses
#' @description List responses with significant associations to a given sample group.
#' @param model NetResponseModel object.
#' @param sample User-specified samples group for which the enrichments are
#' calculated. For instance, an annotation category.
#' @param qth q-value threshold for enrichments
#' @param method Enrichment method.
#' @return Statistics of the significantly associated responses.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @seealso response.enrichment
#' @references See citation("netresponse")
#' @keywords utilities
#' @export
#' @examples # 
list.significant.responses <- function (model, sample, qth = 1, 
			                method = "hypergeometric") {

  # Order responses according to their association with the given sample group
  o <- order.responses(model, sample = sample, method = "hypergeometric")$ordered.responses
  o[which(o$qvalue < qth),]

}




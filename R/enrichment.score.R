enrichment.score <- function(total.samples, subset.samples, annotated.samples, method = "hypergeometric") {
    
    if (method == "hypergeometric") {
        
        N <- length(total.samples)
        
        # number of white balls in the urn
        m <- length(annotated.samples)
        
        # number of black balls in the urn
        n <- N - m
        
        # number of balls drawn from the urn
        k <- length(subset.samples)
        
        # overlap between investigated sample group among response samples (using hard
        # assignments) number of white balls drawn without replacement
        q <- sum(annotated.samples %in% subset.samples)
        
        # hypergeometric enrichment (small p, high enrichment) take 1-p to indicate high
        # enrichment with high score use q-1 since lower.tail = FALSE indicates X > x
        # calculation, but we need X >=x enr <- 1 - phyper(q-1, m, n, k, lower.tail =
        # FALSE, log.p = FALSE)
        pval <- phyper(q - 1, m, n, k, lower.tail = FALSE, log.p = FALSE)
        
        temp <- c(sample.size.total = N, sample.size.subset = k, sample.size.annotated = m, 
            annotated.in.subset = q, fraction.in.data = m/N, fraction.in.subset = q/k, 
            pvalue = pval)
        
        enr <- list(score = 1 - pval, info = temp)
    }
    
    enr
    
}




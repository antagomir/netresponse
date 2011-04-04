ICMg.get.comp.memberships <- function(links, samples) {

  Nnodes <- max(links)
  Nlinks <- dim(links)[1]
  Ncomps <- max(samples$z)
  Nsamples <- dim(samples$z)[1]

  ## Node proportion matrix (average over counts)
  proportions <- matrix(0, Ncomps, Nnodes)

  for (s in 1:Nsamples) {
 
    counts <- matrix(0, Ncomps, Nnodes)
    ## Counts from link components
    for (n in 1:Nlinks) {

      a <- links[n,1]
      b <- links[n,2]
      c <- samples$z[s,n]
      counts[c,a] <- counts[c,a]+1
      counts[c,b] <- counts[c,b]+1
    }

    ## Counts from node components (if given)
    if (length(samples$w)>0) {
      for (m in 1:Nnodes) {
        c <- samples$w[s,m]
        counts[c,m] <- counts[c,m]+1
      }
    }
    
    ## Compute component-wise average from the sample s
    proportions <- proportions + t(t(counts)/apply(counts, 2, sum))
  }
  return(proportions/Nsamples)
}


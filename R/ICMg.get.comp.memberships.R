#
# Copyright (C) 2008-2011 Juuso Parkkinen
# Contact: Juuso Parkkinen <juuso.parkkinen@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

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
    
    ## Compute component-wise average from the samples
    proportions <- proportions + t(t(counts)/apply(counts, 2, sum))
  }
  return(proportions/Nsamples)
}


# Copyright (C) 2010-2012 Leo Lahti
# Contact: Leo Lahti <leo.lahti@iki.fi>
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


#' Reading network files
#' 
#' Function to read network files.
#' 
#' Read in SIF network file, return R graph object in igraph or graphNEL
#' format.
#' 
#' @aliases read.sif
#' @usage read.sif(sif.file, format = "graphNEL", directed = FALSE, header =
#' TRUE, sep = "\t", ...)
#' @param sif.file Name of network file in SIF format.
#' @param format Output format: igraph0 or graphNEL
#' @param directed Logical. Directed/undirected graph. Not used in the current
#' model.
#' @param header Logical. Indicate whether the SIF file has header or not.
#' @param sep Field separator.
#' @param ... Further optional arguments to be passed for file reading.
#' @return R graph object in igraph0 or graphNEL format.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords utilities
#' @export
#' @examples
#' 
#' #net <- read.sif("network.sif")
#' 
read.sif <- function (sif.file, format = "graphNEL", directed = FALSE, header = TRUE, sep = "\t", ...) 
{

    net <- read.csv(file = sif.file, sep = sep, colClasses = "character", header = header, ...)
        
    # Assume form: node1 linktype node2 side.info..
    if ( ncol(net) > 2 ) { 

      # remove NA nodes 
      nas <- apply(net, 1, function (x) {any(is.na(x[c(1,3)]))})
      if (any(nas)) {
        net <- net[!nas, ]
        warning("NAs removed from network node list, ", sum(nas), " edges removed.")
      }
      
      net <- graph.edgelist(as.matrix(net[, -2]), directed = directed)

    } else if ( ncol(net) == 2 ) { # assume form: node1 node2

      # remove NA nodes 
      nas <- apply(net, 1, function (x) {any(is.na(x))})
      if (any(nas)) {
        net <- net[!nas, ]
        warning("NAs removed from network node list, ", sum(nas), " edges removed.")
      }
      
      net <- graph.edgelist(cbind(net[,1],net[,2]), directed = directed)
    }

    if (format == "graphNEL") { net <- igraph.to.graphNEL(net) }
    # if (format == "igraph") { net <- igraph.from.graphNEL(igraph.to.graphNEL(net)) }

    net
}

#library(graph)
#g <- read.graph(file, format = c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "graphdb", "gml"))
# FIXME: check also these; then convert graph to graphNEL



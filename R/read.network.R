read.sif <- 
function (sif.file, format = "graphNEL", directed = FALSE) 
{

    net <- read.csv(file = sif.file, sep = "\t", colClasses = "character")
        
    # assume form: node1 linktype node2 side.info..
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

    net
}

#library(graph)
#g <- read.graph(file, format = c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "graphdb", "gml"))
# FIXME: check also these; then convert graph to graphNEL



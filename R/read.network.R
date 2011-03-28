read.sif <- function (sif.file, format = "graphNEL", directed = FALSE) 
{
  net <- read.table(file = sif.file, colClasses = "character")
  net <- graph.edgelist(as.matrix(net[, -2]), directed = directed)

  if (format == "graphNEL") {
    net <- igraph.to.graphNEL(net)
  }

  net

}

#library(graph)
#g <- read.graph(file, format = c("edgelist", "pajek", "ncol", "lgl", "graphml", "dimacs", "graphdb", "gml"))
# FIXME: check also these; then convert graph to graphNEL



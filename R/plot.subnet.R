plot.subnet <-
function (x, subnet.id, network, plot.names = TRUE, ...) {

  require(Rgraphviz)
  require(igraph)

  if (is.numeric(subnet.id)) {
    subnet.id <- paste("Subnet", subnet.id, sep = "-")
    warning("subnet.id given as numeric; converting to character: ", "Subnet-", subnet.id, sep="")
  }

  subnet.nodes <- get.subnets(x)[[subnet.id]]
  mynet <- network[subnet.nodes, subnet.nodes]

  tmp <- plot.response(x = NULL, mynet, mybreaks = NULL, mypalette = NULL, colors = FALSE, maintext = subnet.id)

  mynet

}

 
 #message("convert to matrix graph format")
 #myg <- new("graphAM", mynet, "undirected")
 #myg2 <-as(myg, "graphNEL") 


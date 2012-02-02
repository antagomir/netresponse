plot.response <-
function (x, mynet, mybreaks, mypalette, plot.names = TRUE, colors = TRUE, plot.type = "twopi",
           ...) {

   # Add node color for specific nodes
   nAttrs <- list()
   if (colors) {
     bins <- check.bins(x, mybreaks)
     nAttrs$fillcolor <- mypalette(length(mybreaks) + 1)[bins]
   } else {
     nAttrs$fillcolor <- rep("white", nrow(mynet))
   }
   names(nAttrs$fillcolor) <- rownames(mynet)

   # add node names for all nodes
   if (plot.names) {
     nodenames <- rownames(mynet)
   } else {
     nodenames <- rep("", nrow(mynet))
   }

   nAttrs$label <- nodenames
   names(nAttrs$label) <- rownames(mynet)

   myg <- as(new("graphAM", mynet, "undirected"), "graphNEL")
   
   plot(myg, y = plot.type, nodeAttrs = nAttrs, ...)

 }


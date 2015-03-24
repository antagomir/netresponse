#' @import dmt
#' @import graph
#' @import parallel
#' @import RColorBrewer
#' @import Rgraphviz
#' @importFrom igraph igraph.to.graphNEL
#' @importFrom igraph igraph.from.graphNEL
#' @importFrom igraph subgraph
#' @importFrom igraph graph.data.frame
#' @importFrom igraph graph.edgelist
#' @import qvalue
#' @import mclust
#' @import methods
#' @importFrom plyr ddply
#' @importFrom minet build.mim
#' @importFrom reshape melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 theme_set
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 opts



.onAttach <- function(lib, pkg)
{
   packageStartupMessage('\nnetresponse (C) 2008-2015 Leo Lahti et al.\n\nhttps://github.com/antagomir/netresponse')
}

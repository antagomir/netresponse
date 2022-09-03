#' @import methods
#' @importFrom BiocStyle pkg_ver
#' @import graph
#' @importFrom graphics axis
#' @importFrom graphics barplot
#' @importFrom graphics hist
#' @importFrom graphics image
#' @importFrom graphics legend
#' @importFrom graphics par
#' @importFrom graphics points
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices gray
#' @importFrom grDevices palette
#' @importFrom grDevices rainbow
#' @import parallel
#' @import RColorBrewer
#' @import Rgraphviz
#' @importFrom igraph igraph.to.graphNEL
#' @importFrom igraph igraph.from.graphNEL
#' @importFrom igraph subgraph
#' @importFrom igraph graph.data.frame
#' @importFrom igraph graph.edgelist
#' @import mclust
#' @importFrom plyr ddply
#' @importFrom minet build.mim
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 qplot
#' @importFrom ggplot2 theme_set
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 facet_wrap
#' @importFrom stats cor
#' @importFrom stats dnorm
#' @importFrom stats na.omit
#' @importFrom stats rgamma
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats p.adjust
#' @importFrom stats phyper
#' @importFrom stats princomp
#' @importFrom stats prcomp
#' @importFrom stats qchisq
#' @importFrom stats sd
#' @importFrom stats t.test
#' @importFrom utils read.csv

.onAttach <- function(lib, pkg) {
    packageStartupMessage("\nnetresponse (C) 2008-2022 Leo Lahti et al.\n\nhttps://github.com/antagomir/netresponse")
}

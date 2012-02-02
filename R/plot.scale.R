plot.scale <-
function (x, y, m = NULL, cex.axis = 1.5, label.step = 2, interval = .1, two.sided = TRUE, label.start = 1, Nlab = 3, ...) {

  require(Rgraphviz)
  #require(igraph)
  
  if (two.sided) {
    
    if (length(m) > 0) {
      x <- set.breaks(m, interval)
    } else {
      mm <- max(x[-c(1, length(x))])
      m <- mm - interval/2
    }
  
    image(t(as.matrix(seq(-mm, mm, length = 100))), col = y(length(x) - 1),
          xaxt = 'n', yaxt = 'n', zlim = range(x), breaks = x)
    
    ndigits <- nchar(unlist(strsplit(as.character(mm), "\\."))[[2]])
    digit.step <- 10^(-ndigits)
    labs <- seq(-mm, mm, by = digit.step)
    start.position <- match(-label.start, round(labs, ndigits))
    end.position <- match(label.start, round(labs, ndigits))
    inds <- seq(start.position, end.position,length = Nlab)
      
    axis(2, at = inds/length(labs),
         labels = labs[inds], cex.axis = cex.axis, las = 2)
  }

  if (!two.sided) {

    mm <- max(x) + 1e6 # infty
    m <- max(x)
 
    labs <- seq(0, m, label.step)
    inds <- sapply(labs,function(lab){min(which(lab<=x))})
  
    image(t(as.matrix(seq(0, m, length = 100))), col = y(length(x) - 1),
          xaxt='n', yaxt='n', zlim = range(x), breaks = x)
    
    int <- 1/(length(x)-1)
    axis(2, at = seq(0, 1, by = int)[inds], labels = labs,
         cex.axis = cex.axis, las = 2)
  }
  
}


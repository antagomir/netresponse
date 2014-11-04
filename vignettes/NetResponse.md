# netresponse - probabilistic tools for functional network analysis

For bug reports and maintainer contact details, see the [README](../README.md) file

## Background 

Condition-specific network activation is characteristic for cellular
systems and other real-world interaction networks. If measurements of
network states are available across a versatile set of conditions or
time points, it becomes possible to construct a global view of network
activation patterns. Different parts of the network respond to
different conditions, and in different ways. Systematic, data-driven
identification of these responses will help to obtain a holistic view
of network activity
[[1](http://bioinformatics.oxfordjournals.org/content/26/21/2713.short)-[2](http://lib.tkk.fi/Diss/2010/isbn9789526033686/)]. This
package provides robust probabilistic algorithms for functional
network analysis
[[1](http://bioinformatics.oxfordjournals.org/content/26/21/2713.short),
[3](http://www.biomedcentral.com/1752-0509/4/4)].

The methods are based on nonparametric probabilistic modeling and
variational learning, and provide general exploratory tools to
investigate the structure ([ICMg](http://www.biomedcentral.com/1752-0509/4/4)) and
context-specific behavior ([NetResponse](http://bioinformatics.oxfordjournals.org/content/26/21/2713.short)) of
interaction networks.  ICMg is used to identify community structure in
interaction networks; NetResponse detects and characterizes
subnetworks that exhibit context-specific activation patterns across
versatile collections of functional measurements, such as gene
expression data. The implementations are partially based on the
agglomerative independent variable group analysis ([AIVGA](http://www.sciencedirect.com/science/article/pii/S0925231208000659))
and variational Dirichlet process Gaussian mixture models
([Kurihara et al. 2007](http://machinelearning.wustl.edu/mlpapers/paper_files/NIPS2006_248.pdf)). The tools are particularly useful for global
exploratory analysis of genome-wide interaction networks and versatile
collections of gene expression data.


## Usage examples

Examples on running NetResponse algorithm and visualizing the
results. The algorithm combines network and functional information to
detect coherent subnetworks that reveal distinct activation modes
across conditions. Kindly cite [this
article](http://bioinformatics.oxfordjournals.org/content/26/21/2713.short).


```r
library(netresponse)

# Generate simulated data
res <- generate.toydata(Dim = 3, Nc = 3, Ns = 200, sd0 = 3, rgam.shape = 1, rgam.scale = 1, rseed = 123456)

D <- res$data
component.means <- res$means
component.sds   <- res$sds
sample2comp     <- res$sample2comp

# Fit NetResponse model
# Various network formats are supported, see help(detect.responses) for
# details. With large data sets, consider the 'speedup' option.
set.seed(4243)
res <- detect.responses(D, mixture.method = "vdp", pca.basis = TRUE)

# List subnets (each is a list of nodes)
subnet.id <- names(get.subnets(res))[[1]]
```

### PCA visualization


```r
library(ggplot2)
vis <- plot.responses(res, subnet.id, plot.mode = "pca")
```

![plot of chunk NetResponse2](figure/NetResponse2.png) 


```r
# Modify the resulting ggplot2 object to enhance visualization
p <- vis$p # Pick the ggplot2 object from results
p <- p + geom_point(size = 3) # Modify point size
print(p) # Plot
```

![plot of chunk NetResponse2b](figure/NetResponse2b.png) 


### Network visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "network")
```

![plot of chunk NetResponse3](figure/NetResponse3.png) 

### Heatmap visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "heatmap")
```

```
## [1] "Here"
```

![plot of chunk NetResponse4](figure/NetResponse4.png) 

### Boxplot visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "boxplot.data")
```

![plot of chunk NetResponse5](figure/NetResponse5.png) 

See also mode = "response.barplot" 


### Color scale


```r
plot.scale(vis$breaks, vis$palette, two.sided = TRUE)
```

![plot of chunk NetResponse7](figure/NetResponse7.png) 


### Cluster assignments

The sample-response assignments from the mixture model are soft
ie. defined as continuous probabilities. Retrieve the hard clustering
ie. list of samples for each response, response for each sample, based
the highest probability:


```r
subnet.id <- 'Subnet-1'
sample.probs <- response2sample(res, subnet.id)
response.probs <- sample2response(res, subnet.id)
```

Retrieve model parameters for a given subnetwork (Gaussian mixture
means, covariance diagonal, and component weights):


```r
params <- get.model.parameters(res, subnet.id) 
names(params)
```

```
## [1] "mu"          "sd"          "w"           "free.energy" "Nparams"    
## [6] "qofz"        "nodes"
```



## Nonparametric Gaussian mixture models

The package provides additional tools for nonparametric Gaussian
mixture modeling based on variational Dirichlet process mixture models
and implementations by [Kurihara et al.](http://machinelearning.wustl.edu/mlpapers/paper_files/NIPS2006_248.pdf) and [Honkela et al.](http://www.sciencedirect.com/science/article/pii/S0925231208000659). See the
example in help(vdp.mixt).

## Interaction Component Model for Gene Modules

Interaction Component Model ([ICMg](http://www.biomedcentral.com/1752-0509/4/4)) can be used to find functional gene
modules from either protein interaction data or from combinations of
protein interaction and gene expression data. Run ICMg and cluster the
nodes:


```r
library(netresponse)
data(osmo)
res <- ICMg.combined.sampler(osmo$ppi, osmo$exp, C=10)
```

```
## Sampling ICMg2...
## nodes: 10250 links: 1711 observations: 133 components: 10 alpha: 10 beta: 0.01 
## Sampling 1000 iterationcs
## Burnin iterations: 800 
## I: 0
## n(z): 1043 996 990 1028 1013 1050 995 1068 1023 1044 
## m(z): 166 185 185 149 167 184 160 169 179 167 
## I: 100 
## convL: -0.2279 n(z): 368 759 478 482 599 1071 3457 595 1883 558 
## convN: -0.002489 m(z): 69 169 114 176 103 234 340 168 218 120 
## I: 200 
## convL: -0.2226 n(z): 359 699 425 529 586 1114 3550 659 1758 571 
## convN: -0.001637 m(z): 69 167 115 176 100 235 341 169 220 119 
## I: 300 
## convL: -0.2296 n(z): 356 751 452 506 611 1051 3486 639 1846 552 
## convN: -0.004019 m(z): 68 169 114 176 101 236 338 169 221 119 
## I: 400 
## convL: -0.2226 n(z): 329 747 443 528 584 1110 3337 588 1982 602 
## convN: -0.006692 m(z): 69 168 115 174 101 235 339 170 221 119 
## I: 500 
## convL: -0.2041 n(z): 323 774 424 499 617 1171 3427 603 1883 529 
## convN: -0.004346 m(z): 69 164 116 175 101 235 343 169 220 119 
## I: 600 
## convL: -0.2088 n(z): 356 735 427 524 627 1088 3416 646 1837 594 
## convN: -0.001701 m(z): 68 164 115 176 100 235 345 168 220 120 
## I: 700 
## convL: -0.2111 n(z): 342 779 422 482 638 1172 3365 614 1899 537 
## convN: -0.006655 m(z): 68 167 114 175 101 235 342 170 220 119 
## I: 800 
## convL: -0.222 n(z): 334 792 407 498 589 1129 3357 630 1955 559 
## convN: -0.004293 m(z): 67 166 121 177 96 235 341 168 219 121 
## Sample iterations: 200 
## I: 810 
## convL: -0.2215 n(z): 341 764 440 459 592 1080 3447 669 1930 528 
## convN: -0.002232 m(z): 68 167 120 177 98 232 336 171 222 120 
## I: 820 
## convL: -0.2168 n(z): 350 753 421 468 630 1054 3477 661 1871 565 
## convN: -0.005799 m(z): 68 168 121 176 99 232 335 171 222 119 
## I: 830 
## convL: -0.2231 n(z): 339 760 442 504 574 1093 3434 606 1947 551 
## convN: -0.001073 m(z): 67 168 120 178 96 234 337 170 221 120 
## I: 840 
## convL: -0.2224 n(z): 345 763 419 510 575 1090 3406 603 2006 533 
## convN: -0.005693 m(z): 68 166 121 177 98 232 339 170 221 119 
## I: 850 
## convL: -0.2215 n(z): 366 765 418 475 599 1132 3461 604 1891 539 
## convN: -0.002615 m(z): 67 174 121 179 98 235 342 170 206 119 
## I: 860 
## convL: -0.2089 n(z): 369 708 429 517 622 1113 3370 593 1944 585 
## convN: -0.004933 m(z): 67 175 121 177 97 236 342 171 206 119 
## I: 870 
## convL: -0.2212 n(z): 378 733 437 522 591 1203 3366 597 1843 580 
## convN: -0.002067 m(z): 67 173 121 177 97 236 346 171 204 119 
## I: 880 
## convL: -0.2136 n(z): 356 755 402 509 606 1164 3438 584 1896 540 
## convN: -0.00464 m(z): 67 174 120 177 97 236 345 171 205 119 
## I: 890 
## convL: -0.2312 n(z): 356 839 366 503 597 1201 3357 565 1933 533 
## convN: -0.002128 m(z): 67 173 121 177 98 232 344 171 208 120 
## I: 900 
## convL: -0.2161 n(z): 346 839 357 557 613 1138 3491 527 1806 576 
## convN: -0.006676 m(z): 67 174 121 179 98 235 344 170 205 118 
## I: 910 
## convL: -0.2204 n(z): 351 825 380 527 584 1223 3407 579 1842 532 
## convN: -0.003132 m(z): 67 173 121 177 99 235 344 171 205 119 
## I: 920 
## convL: -0.2073 n(z): 348 813 346 548 587 1141 3519 624 1778 546 
## convN: -0.004245 m(z): 67 174 121 175 99 231 344 172 209 119 
## I: 930 
## convL: -0.2111 n(z): 339 838 373 509 586 1208 3464 563 1834 536 
## convN: -0.002986 m(z): 67 173 121 177 98 234 343 171 207 120 
## I: 940 
## convL: -0.219 n(z): 351 822 386 567 600 1154 3454 554 1837 525 
## convN: -0.002636 m(z): 67 174 120 178 97 236 346 170 204 119 
## I: 950 
## convL: -0.207 n(z): 363 786 368 501 585 1294 3404 551 1848 550 
## convN: -0.005655 m(z): 67 174 121 177 96 237 344 171 204 120 
## I: 960 
## convL: -0.218 n(z): 326 727 422 556 568 1194 3505 579 1831 542 
## convN: -0.001847 m(z): 67 173 121 177 98 235 344 171 205 120 
## I: 970 
## convL: -0.2078 n(z): 333 762 392 515 580 1264 3488 557 1806 553 
## convN: -0.003445 m(z): 67 173 121 177 97 234 346 171 205 120 
## I: 980 
## convL: -0.2149 n(z): 340 771 364 542 588 1145 3510 615 1800 575 
## convN: -0.002504 m(z): 67 173 121 175 100 234 345 172 204 120 
## I: 990 
## convL: -0.2176 n(z): 342 781 397 522 613 1228 3445 569 1763 590 
## convN: -0.001913 m(z): 68 175 121 177 97 236 342 171 205 119 
## I: 1000 
## convL: -0.2206 n(z): 330 773 380 592 576 1202 3470 576 1805 546 
## convN: -0.003711 m(z): 67 174 121 176 99 236 344 171 204 119 
## DONE
```

```r
res$comp.memb <- ICMg.get.comp.memberships(osmo$ppi, res)
res$clustering <- apply(res$comp.memb, 2, which.max)
```


### Citing NetResponse

Please cite [Lahti et al. (2010)](http://bioinformatics.oxfordjournals.org/content/26/21/2713) with the package. When using the ICMg algorithms, additionally cite [Parkkinen et al. (2010)](http://www.biomedcentral.com/1752-0509/4/4).


### Version information

This document was written using:


```r
sessionInfo()
```

```
## R version 3.1.2 (2014-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] ggplot2_1.0.0       netresponse_1.17.12 reshape_0.8.5      
## [4] mclust_4.4          minet_3.20.1        Rgraphviz_2.8.1    
## [7] graph_1.42.0        knitr_1.6          
## 
## loaded via a namespace (and not attached):
##  [1] BiocGenerics_0.10.0 colorspace_1.2-4    digest_0.6.4       
##  [4] dmt_0.8.20          evaluate_0.5.5      formatR_1.0        
##  [7] gtable_0.1.2        igraph_0.7.1        labeling_0.3       
## [10] lattice_0.20-29     MASS_7.3-34         Matrix_1.1-4       
## [13] munsell_0.4.2       mvtnorm_1.0-0       parallel_3.1.2     
## [16] plyr_1.8.1          proto_0.3-10        qvalue_1.38.0      
## [19] RColorBrewer_1.0-5  Rcpp_0.11.2         reshape2_1.4       
## [22] scales_0.2.4        stats4_3.1.2        stringr_0.6.2      
## [25] tcltk_3.1.2         tools_3.1.2
```

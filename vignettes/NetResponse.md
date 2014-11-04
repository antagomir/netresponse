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

![plot of chunk NetResponse2](figure/NetResponse21.png) 

```r
# Modify the resulting ggplot2 object to enhance visualization
vis$p + geom_point(size = 3)
```

![plot of chunk NetResponse2](figure/NetResponse22.png) 

### Network visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "network")
```

```
## Error: could not find function "check.bins"
```

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

```
## Error: could not find function "get.dat"
```

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
## convL: -0.2397 n(z): 282 603 828 691 475 1038 3542 727 1572 492 
## convN: -0.005183 m(z): 85 105 183 164 97 199 382 183 185 128 
## I: 200 
## convL: -0.2372 n(z): 337 631 822 632 413 999 3410 804 1743 459 
## convN: -0.002851 m(z): 85 104 186 163 97 195 382 183 187 129 
## I: 300 
## convL: -0.2295 n(z): 345 591 850 668 385 966 3655 678 1653 459 
## convN: -0.0009674 m(z): 86 105 183 164 98 195 385 180 186 129 
## I: 400 
## convL: -0.216 n(z): 367 551 846 602 399 963 3663 614 1775 470 
## convN: -0.001825 m(z): 86 105 183 164 98 197 382 180 187 129 
## I: 500 
## convL: -0.2393 n(z): 350 656 859 569 436 1011 3574 524 1779 492 
## convN: -0.005055 m(z): 85 107 182 167 97 194 385 182 186 126 
## I: 600 
## convL: -0.2226 n(z): 353 680 840 501 444 961 3590 496 1886 499 
## convN: -0.001065 m(z): 86 107 183 165 97 197 380 180 188 128 
## I: 700 
## convL: -0.2274 n(z): 341 585 824 605 465 1067 3583 507 1812 461 
## convN: -0.001451 m(z): 86 107 185 167 98 195 381 180 186 126 
## I: 800 
## convL: -0.2251 n(z): 348 623 828 585 430 1063 3580 584 1717 492 
## convN: -0.003788 m(z): 86 105 184 165 100 203 379 180 184 125 
## Sample iterations: 200 
## I: 810 
## convL: -0.2233 n(z): 327 599 819 605 373 1072 3605 515 1849 486 
## convN: -0.005114 m(z): 69 108 189 184 102 204 380 177 181 117 
## I: 820 
## convL: -0.2238 n(z): 371 565 784 637 368 1112 3598 543 1784 488 
## convN: -0.002216 m(z): 71 107 189 184 104 203 379 176 182 116 
## I: 830 
## convL: -0.2237 n(z): 352 577 863 580 364 1072 3568 550 1825 499 
## convN: -0.002508 m(z): 70 108 190 183 100 206 381 177 180 116 
## I: 840 
## convL: -0.2262 n(z): 340 595 846 599 359 1107 3503 551 1859 491 
## convN: -0.002692 m(z): 71 107 190 184 99 206 380 177 182 115 
## I: 850 
## convL: -0.2184 n(z): 344 603 842 594 407 1076 3572 557 1782 473 
## convN: -0.001712 m(z): 71 107 190 183 100 205 380 176 183 116 
## I: 860 
## convL: -0.2177 n(z): 326 568 864 570 384 993 3582 568 1868 527 
## convN: -0.002416 m(z): 70 108 190 183 100 206 381 177 180 116 
## I: 870 
## convL: -0.2205 n(z): 330 607 838 552 411 1059 3596 576 1787 494 
## convN: -0.001688 m(z): 71 108 190 183 100 206 381 177 180 115 
## I: 880 
## convL: -0.2185 n(z): 359 644 833 542 410 1087 3591 528 1784 472 
## convN: -0.0009068 m(z): 71 108 190 183 100 206 381 177 180 115 
## I: 890 
## convL: -0.2203 n(z): 355 590 840 535 408 1065 3555 534 1916 452 
## convN: -0.006994 m(z): 71 107 191 183 100 206 377 177 184 115 
## I: 900 
## convL: -0.2273 n(z): 348 588 884 544 401 1100 3532 550 1813 490 
## convN: -0.001342 m(z): 71 107 190 183 100 206 380 176 182 116 
## I: 910 
## convL: -0.2259 n(z): 377 573 852 609 410 1014 3576 508 1836 495 
## convN: -0.004068 m(z): 70 107 190 183 100 207 379 176 182 117 
## I: 920 
## convL: -0.216 n(z): 337 543 864 605 402 1039 3556 538 1855 511 
## convN: -0.007247 m(z): 70 107 189 183 102 206 381 175 181 117 
## I: 930 
## convL: -0.2274 n(z): 318 547 854 678 402 1128 3539 494 1805 485 
## convN: -0.003043 m(z): 69 108 190 183 100 204 380 176 183 118 
## I: 940 
## convL: -0.2155 n(z): 351 536 767 658 417 1113 3545 507 1875 481 
## convN: -0.006344 m(z): 70 108 190 183 100 204 379 176 184 117 
## I: 950 
## convL: -0.2028 n(z): 368 536 815 675 394 1074 3597 533 1797 461 
## convN: -0.003304 m(z): 71 108 190 184 100 207 381 174 179 117 
## I: 960 
## convL: -0.2287 n(z): 316 578 892 587 401 1127 3557 533 1809 450 
## convN: -0.001386 m(z): 70 107 190 183 100 207 381 177 180 116 
## I: 970 
## convL: -0.2222 n(z): 316 587 870 610 415 1050 3562 531 1813 496 
## convN: -0.00115 m(z): 71 107 190 183 100 206 381 176 181 116 
## I: 980 
## convL: -0.2191 n(z): 355 588 844 637 390 1060 3601 488 1815 472 
## convN: -0.002459 m(z): 71 109 190 183 100 202 381 176 183 116 
## I: 990 
## convL: -0.2264 n(z): 309 582 882 605 428 1058 3528 509 1861 488 
## convN: -0.003876 m(z): 70 108 186 183 92 207 382 177 188 118 
## I: 1000 
## convL: -0.2158 n(z): 321 609 893 617 382 1030 3573 554 1814 457 
## convN: -0.001571 m(z): 71 109 188 183 91 208 381 177 186 117 
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

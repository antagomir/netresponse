# netresponse - probabilistic tools for functional network analysis

For maintainer contact details, see the [README](../README.md) file

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
vis <- plot.responses(res, subnet.id, plot.mode = "pca")
vis$p + geom_point(size = 3)
```

![plot of chunk NetResponse2](figure/NetResponse21.png) ![plot of chunk NetResponse2](figure/NetResponse22.png) 

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
## Error: could not find function "plotMatrix.2way"
```

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
## n(z): 1044 995 991 1028 1013 1050 995 1068 1023 1043 
## m(z): 166 185 185 149 167 184 160 168 179 168 
## I: 100 
## convL: -0.2275 n(z): 507 651 806 419 1045 595 1725 596 473 3433 
## convN: -0.001644 m(z): 127 165 200 168 204 121 185 109 100 332 
## I: 200 
## convL: -0.2333 n(z): 462 675 722 547 956 582 1865 572 434 3435 
## convN: -0.003056 m(z): 135 168 190 172 200 120 184 113 97 332 
## I: 300 
## convL: -0.2165 n(z): 467 660 716 499 1060 581 1810 555 442 3460 
## convN: -0.001647 m(z): 135 167 191 170 198 119 186 114 97 334 
## I: 400 
## convL: -0.2298 n(z): 448 681 675 475 1053 598 1967 526 420 3407 
## convN: -0.002638 m(z): 135 167 191 170 200 120 186 114 96 332 
## I: 500 
## convL: -0.2251 n(z): 476 649 692 499 1019 583 1885 569 407 3471 
## convN: -0.004271 m(z): 135 168 191 170 201 120 181 115 96 334 
## I: 600 
## convL: -0.2142 n(z): 495 684 742 526 905 594 1909 526 408 3461 
## convN: -0.002076 m(z): 135 168 191 170 199 120 185 114 97 332 
## I: 700 
## convL: -0.2268 n(z): 525 726 662 464 978 622 1887 515 432 3439 
## convN: -0.002185 m(z): 136 168 190 169 200 120 184 114 98 332 
## I: 800 
## convL: -0.2291 n(z): 496 754 700 470 941 602 1913 496 429 3449 
## convN: -0.00243 m(z): 136 167 191 169 199 120 188 114 96 331 
## Sample iterations: 200 
## I: 810 
## convL: -0.2351 n(z): 493 721 705 496 950 570 1934 523 404 3454 
## convN: -0.001951 m(z): 137 168 190 169 200 121 183 114 97 332 
## I: 820 
## convL: -0.2305 n(z): 473 675 677 524 1018 565 1893 558 434 3433 
## convN: -0.005067 m(z): 136 167 191 169 199 120 184 113 97 335 
## I: 830 
## convL: -0.2276 n(z): 469 705 718 471 970 618 1871 538 412 3478 
## convN: -0.003448 m(z): 136 168 190 169 200 121 183 114 97 333 
## I: 840 
## convL: -0.2266 n(z): 468 678 719 487 979 601 1939 536 405 3438 
## convN: -0.001018 m(z): 136 167 191 169 200 120 183 115 97 333 
## I: 850 
## convL: -0.2147 n(z): 520 621 701 439 1028 618 1921 553 406 3443 
## convN: -0.004278 m(z): 137 168 191 168 200 120 184 114 98 331 
## I: 860 
## convL: -0.2186 n(z): 473 673 693 475 950 621 1975 521 399 3470 
## convN: -0.001615 m(z): 136 168 191 169 199 120 185 114 97 332 
## I: 870 
## convL: -0.2164 n(z): 465 700 682 536 981 640 1817 519 406 3504 
## convN: -0.004295 m(z): 136 168 190 168 198 121 185 113 98 334 
## I: 880 
## convL: -0.2111 n(z): 491 685 669 501 970 619 1817 508 435 3555 
## convN: -0.005394 m(z): 136 168 190 170 199 119 183 115 97 334 
## I: 890 
## convL: -0.2047 n(z): 515 712 680 513 978 603 1768 525 467 3489 
## convN: -0.003229 m(z): 136 167 191 169 199 120 186 114 97 332 
## I: 900 
## convL: -0.2201 n(z): 520 674 694 497 983 588 1928 498 417 3451 
## convN: -0.003665 m(z): 136 168 190 170 200 119 184 114 97 333 
## I: 910 
## convL: -0.2117 n(z): 506 720 645 490 944 599 1988 507 432 3419 
## convN: -0.001064 m(z): 136 169 190 169 200 120 182 114 98 333 
## I: 920 
## convL: -0.2186 n(z): 479 704 673 492 966 551 1956 519 428 3482 
## convN: -0.002769 m(z): 136 168 190 170 199 120 185 114 96 333 
## I: 930 
## convL: -0.2095 n(z): 489 716 704 526 973 600 1828 525 418 3471 
## convN: -0.002302 m(z): 136 169 191 169 200 119 185 114 97 331 
## I: 940 
## convL: -0.2325 n(z): 467 672 721 508 980 590 1842 610 401 3459 
## convN: -0.001974 m(z): 136 167 190 169 200 121 184 114 97 333 
## I: 950 
## convL: -0.2252 n(z): 499 692 641 555 1015 629 1828 579 405 3407 
## convN: -0.002789 m(z): 136 168 190 169 200 121 185 114 97 331 
## I: 960 
## convL: -0.2301 n(z): 478 693 704 535 1029 624 1759 573 387 3468 
## convN: -0.002249 m(z): 136 169 190 169 199 120 185 114 98 331 
## I: 970 
## convL: -0.2268 n(z): 504 669 673 504 936 589 1859 552 448 3516 
## convN: -0.001774 m(z): 137 167 190 169 201 120 182 114 98 333 
## I: 980 
## convL: -0.2275 n(z): 509 706 730 511 989 619 1808 554 392 3432 
## convN: -0.008113 m(z): 136 168 190 170 201 119 182 114 96 335 
## I: 990 
## convL: -0.2176 n(z): 533 662 697 499 959 627 1811 562 403 3497 
## convN: -0.0006697 m(z): 136 168 190 170 199 119 184 114 98 333 
## I: 1000 
## convL: -0.2146 n(z): 464 651 698 472 922 604 1927 568 419 3525 
## convN: -0.001428 m(z): 136 169 190 169 198 120 183 114 98 334 
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
## R version 3.1.1 (2014-07-10)
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
## [1] ggplot2_1.0.0      netresponse_1.17.1 reshape_0.8.5     
## [4] mclust_4.4         minet_3.20.1       igraph_0.7.1      
## [7] Rgraphviz_2.8.1    graph_1.42.0       knitr_1.6         
## 
## loaded via a namespace (and not attached):
##  [1] BiocGenerics_0.10.0 colorspace_1.2-4    digest_0.6.4       
##  [4] dmt_0.8.20          evaluate_0.5.5      formatR_1.0        
##  [7] gtable_0.1.2        labeling_0.3        lattice_0.20-29    
## [10] MASS_7.3-34         Matrix_1.1-4        munsell_0.4.2      
## [13] mvtnorm_1.0-0       parallel_3.1.1      plyr_1.8.1         
## [16] proto_0.3-10        qvalue_1.38.0       RColorBrewer_1.0-5 
## [19] Rcpp_0.11.2         reshape2_1.4        scales_0.2.4       
## [22] stats4_3.1.1        stringr_0.6.2       tcltk_3.1.1        
## [25] tools_3.1.1
```

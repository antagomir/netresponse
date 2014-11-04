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
## n(z): 1044 995 991 1028 1013 1050 995 1068 1023 1043 
## m(z): 166 185 185 149 167 184 160 168 179 168 
## I: 100 
## convL: -0.2675 n(z): 629 696 802 568 899 481 1498 877 479 3321 
## convN: -0.004126 m(z): 137 167 191 169 201 121 178 116 98 333 
## I: 200 
## convL: -0.2413 n(z): 562 729 793 629 812 452 1510 857 448 3458 
## convN: -0.004319 m(z): 137 169 192 170 198 120 178 117 97 333 
## I: 300 
## convL: -0.2288 n(z): 589 756 776 643 824 411 1475 849 466 3461 
## convN: -0.002966 m(z): 137 171 192 170 195 120 176 118 97 335 
## I: 400 
## convL: -0.2302 n(z): 596 820 746 688 883 416 1368 822 446 3465 
## convN: -0.001239 m(z): 137 170 193 168 195 120 178 118 98 334 
## I: 500 
## convL: -0.2353 n(z): 585 870 730 636 905 485 1230 927 458 3424 
## convN: -0.005092 m(z): 133 175 190 170 196 122 141 148 96 340 
## I: 600 
## convL: -0.2388 n(z): 598 821 701 620 864 480 1210 1041 470 3445 
## convN: -0.001264 m(z): 134 170 192 169 194 121 140 154 94 343 
## I: 700 
## convL: -0.2262 n(z): 568 802 724 637 874 428 1212 1067 418 3520 
## convN: -0.002102 m(z): 134 169 193 169 194 121 140 154 94 343 
## I: 800 
## convL: -0.2301 n(z): 628 840 731 453 880 446 1246 1096 455 3475 
## convN: -0.00262 m(z): 134 170 191 171 193 121 140 154 95 342 
## Sample iterations: 200 
## I: 810 
## convL: -0.2362 n(z): 588 889 749 521 866 434 1198 1064 453 3488 
## convN: -0.00211 m(z): 134 170 192 170 193 121 140 154 95 342 
## I: 820 
## convL: -0.2305 n(z): 590 811 708 534 873 445 1253 1159 436 3441 
## convN: -0.001315 m(z): 134 168 192 170 196 122 141 154 94 340 
## I: 830 
## convL: -0.2304 n(z): 590 834 703 482 947 446 1277 1043 449 3479 
## convN: -0.002916 m(z): 134 170 193 170 195 120 139 153 94 343 
## I: 840 
## convL: -0.2267 n(z): 572 791 770 482 872 461 1249 1150 448 3455 
## convN: -0.000676 m(z): 134 170 192 170 194 120 141 154 94 342 
## I: 850 
## convL: -0.2328 n(z): 573 778 804 559 855 430 1333 998 428 3492 
## convN: -0.001185 m(z): 134 170 192 170 194 120 140 154 94 343 
## I: 860 
## convL: -0.2267 n(z): 591 804 726 510 881 454 1282 1038 471 3493 
## convN: -0.003309 m(z): 134 169 190 171 194 122 140 154 94 343 
## I: 870 
## convL: -0.2291 n(z): 601 856 704 512 899 444 1292 1067 440 3435 
## convN: -0.002862 m(z): 134 170 192 169 194 121 141 154 94 342 
## I: 880 
## convL: -0.2329 n(z): 620 817 750 553 828 459 1261 1054 469 3439 
## convN: -0.00116 m(z): 134 170 193 170 194 120 140 154 94 342 
## I: 890 
## convL: -0.2367 n(z): 592 817 766 633 817 470 1153 1134 469 3399 
## convN: -0.0007306 m(z): 134 170 190 170 195 122 141 154 94 341 
## I: 900 
## convL: -0.2268 n(z): 605 779 771 645 749 450 1237 1062 471 3481 
## convN: -0.0005595 m(z): 134 170 193 170 193 120 140 154 94 343 
## I: 910 
## convL: -0.2282 n(z): 591 858 752 635 803 427 1251 1042 470 3421 
## convN: -0.00286 m(z): 134 170 191 170 195 121 141 154 94 341 
## I: 920 
## convL: -0.2347 n(z): 558 850 744 635 784 429 1279 1018 495 3458 
## convN: -0.0009814 m(z): 134 170 190 170 194 122 141 154 94 342 
## I: 930 
## convL: -0.2244 n(z): 592 877 723 565 777 470 1272 1121 438 3415 
## convN: -0.001251 m(z): 134 170 191 170 196 121 141 153 94 341 
## I: 940 
## convL: -0.2338 n(z): 603 826 747 644 733 453 1191 1108 463 3482 
## convN: -0.0002326 m(z): 134 171 191 169 194 121 141 154 94 342 
## I: 950 
## convL: -0.2344 n(z): 613 769 718 671 760 468 1170 1078 468 3535 
## convN: -0.001298 m(z): 134 170 190 170 195 121 141 154 94 342 
## I: 960 
## convL: -0.2322 n(z): 595 792 746 679 822 409 1208 1030 481 3488 
## convN: -0.001364 m(z): 134 171 190 170 194 121 141 154 94 342 
## I: 970 
## convL: -0.2357 n(z): 599 807 719 614 835 441 1179 1089 481 3486 
## convN: -0.001572 m(z): 134 170 193 170 194 120 139 154 94 343 
## I: 980 
## convL: -0.2267 n(z): 597 827 747 599 811 442 1206 1116 437 3468 
## convN: -0.00113 m(z): 134 170 192 170 194 120 141 154 94 342 
## I: 990 
## convL: -0.2333 n(z): 577 813 748 623 825 465 1227 1077 471 3424 
## convN: -0.001292 m(z): 135 170 191 169 194 122 141 154 93 342 
## I: 1000 
## convL: -0.222 n(z): 576 818 716 608 789 472 1153 1140 451 3527 
## convN: -0.005577 m(z): 134 168 193 169 194 123 140 153 95 342 
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

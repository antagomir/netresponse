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

```
## Error: could not find function "geom_point"
```

![plot of chunk NetResponse2](figure/NetResponse2.png) 

### Network visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "network")
```

```
## Error: object 'p' not found
```

![plot of chunk NetResponse3](figure/NetResponse3.png) 

### Heatmap visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "heatmap")
```

```
## Error: object 'p' not found
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
## n(z): 1044 995 991 1028 1013 1050 995 1068 1023 1043 
## m(z): 166 185 185 149 167 184 160 168 179 168 
## I: 100 
## convL: -0.2456 n(z): 494 609 819 1008 703 684 1437 407 593 3496 
## convN: -0.0005979 m(z): 128 157 205 187 164 121 149 96 166 338 
## I: 200 
## convL: -0.2371 n(z): 504 632 814 992 740 633 1484 402 552 3497 
## convN: -0.001942 m(z): 129 157 203 187 164 122 149 96 166 338 
## I: 300 
## convL: -0.2281 n(z): 455 630 820 1019 735 547 1423 415 666 3540 
## convN: -0.001252 m(z): 128 157 203 186 161 122 155 97 167 335 
## I: 400 
## convL: -0.2285 n(z): 500 619 738 949 735 553 1543 425 711 3477 
## convN: -0.001022 m(z): 128 157 204 189 161 122 154 97 166 333 
## I: 500 
## convL: -0.2323 n(z): 509 633 749 979 731 614 1494 384 625 3532 
## convN: -0.001 m(z): 127 157 205 188 161 122 155 97 166 333 
## I: 600 
## convL: -0.2354 n(z): 507 639 695 888 730 597 1617 397 648 3532 
## convN: -0.001869 m(z): 129 157 204 189 161 122 154 97 165 333 
## I: 700 
## convL: -0.2287 n(z): 515 665 719 903 692 577 1580 417 741 3441 
## convN: -0.001457 m(z): 127 157 205 188 161 122 155 97 165 334 
## I: 800 
## convL: -0.2293 n(z): 530 671 770 888 748 540 1581 400 609 3513 
## convN: -0.002125 m(z): 129 157 204 188 161 122 154 97 165 334 
## Sample iterations: 200 
## I: 810 
## convL: -0.2282 n(z): 479 578 765 972 739 603 1563 380 606 3565 
## convN: -0.001088 m(z): 128 158 203 188 161 121 155 97 167 333 
## I: 820 
## convL: -0.2214 n(z): 496 582 730 882 759 563 1519 421 696 3602 
## convN: -0.002145 m(z): 127 159 206 187 161 120 155 97 165 334 
## I: 830 
## convL: -0.2192 n(z): 498 561 776 940 732 584 1575 406 664 3514 
## convN: -0.001644 m(z): 127 164 203 187 162 120 154 97 164 333 
## I: 840 
## convL: -0.2346 n(z): 482 599 819 999 684 591 1523 398 681 3474 
## convN: -0.004877 m(z): 129 154 206 201 164 118 159 97 165 318 
## I: 850 
## convL: -0.2282 n(z): 473 522 781 983 746 580 1666 417 576 3506 
## convN: -0.002188 m(z): 131 148 202 197 167 122 163 98 163 320 
## I: 860 
## convL: -0.2437 n(z): 484 526 798 1034 813 574 1513 429 597 3482 
## convN: -0.001439 m(z): 131 149 206 212 165 118 156 96 162 316 
## I: 870 
## convL: -0.2337 n(z): 483 546 698 1039 858 533 1647 407 563 3476 
## convN: -0.001771 m(z): 126 149 207 214 165 119 156 96 163 316 
## I: 880 
## convL: -0.2253 n(z): 497 560 710 973 826 550 1581 419 613 3521 
## convN: -0.002001 m(z): 125 149 208 214 165 119 154 97 164 316 
## I: 890 
## convL: -0.2281 n(z): 480 601 747 1003 910 536 1588 418 514 3453 
## convN: -0.0006936 m(z): 126 149 206 211 166 119 156 97 165 316 
## I: 900 
## convL: -0.2301 n(z): 467 554 787 1049 805 624 1465 400 557 3542 
## convN: -0.003258 m(z): 126 149 208 213 165 119 157 95 163 316 
## I: 910 
## convL: -0.2397 n(z): 477 536 851 1060 921 601 1396 374 520 3514 
## convN: -0.002268 m(z): 127 149 206 213 165 119 157 95 164 316 
## I: 920 
## convL: -0.2394 n(z): 479 556 796 1067 944 547 1445 387 513 3516 
## convN: -0.00109 m(z): 126 150 206 212 165 119 156 97 165 315 
## I: 930 
## convL: -0.2321 n(z): 493 546 762 1029 995 569 1399 396 515 3546 
## convN: -0.004516 m(z): 126 150 205 212 164 120 156 97 165 316 
## I: 940 
## convL: -0.2331 n(z): 469 607 710 1027 916 553 1486 407 532 3543 
## convN: -0.002793 m(z): 125 150 208 215 164 119 154 97 163 316 
## I: 950 
## convL: -0.2303 n(z): 468 537 762 999 972 554 1441 386 619 3512 
## convN: -0.01016 m(z): 126 149 208 211 166 118 155 96 165 317 
## I: 960 
## convL: -0.2406 n(z): 461 496 679 1052 953 588 1475 399 644 3503 
## convN: -0.001481 m(z): 125 149 208 214 165 119 155 97 163 316 
## I: 970 
## convL: -0.2368 n(z): 462 540 737 1090 936 576 1383 401 621 3504 
## convN: -0.001676 m(z): 125 150 207 211 166 119 157 96 165 315 
## I: 980 
## convL: -0.2453 n(z): 465 593 773 991 1009 557 1506 368 567 3421 
## convN: -0.001355 m(z): 126 151 206 212 165 119 156 97 165 314 
## I: 990 
## convL: -0.2234 n(z): 465 538 806 937 988 606 1427 385 535 3563 
## convN: -0.002933 m(z): 125 150 207 213 165 119 155 97 165 315 
## I: 1000 
## convL: -0.2227 n(z): 472 574 734 939 932 610 1533 391 544 3521 
## convN: -0.001361 m(z): 126 150 206 211 166 119 156 97 165 315 
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
## [1] netresponse_1.17.11 reshape_0.8.5       mclust_4.4         
## [4] minet_3.20.1        igraph_0.7.1        Rgraphviz_2.8.1    
## [7] graph_1.42.0        knitr_1.6          
## 
## loaded via a namespace (and not attached):
##  [1] BiocGenerics_0.10.0 colorspace_1.2-4    digest_0.6.4       
##  [4] dmt_0.8.20          evaluate_0.5.5      formatR_1.0        
##  [7] ggplot2_1.0.0       gtable_0.1.2        labeling_0.3       
## [10] lattice_0.20-29     MASS_7.3-34         Matrix_1.1-4       
## [13] munsell_0.4.2       mvtnorm_1.0-0       parallel_3.1.1     
## [16] plyr_1.8.1          proto_0.3-10        qvalue_1.38.0      
## [19] RColorBrewer_1.0-5  Rcpp_0.11.2         reshape2_1.4       
## [22] scales_0.2.4        stats4_3.1.1        stringr_0.6.2      
## [25] tcltk_3.1.1         tools_3.1.1
```

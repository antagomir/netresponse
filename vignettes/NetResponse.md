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
vis$p + geom_point(size = 3)
```

![plot of chunk NetResponse2](figure/NetResponse21.png) ![plot of chunk NetResponse2](figure/NetResponse22.png) 

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
## n(z): 1043 996 990 1028 1013 1050 995 1068 1023 1044 
## m(z): 166 185 185 149 167 184 160 169 179 167 
## I: 100 
## convL: -0.2242 n(z): 219 513 631 689 496 1172 3552 704 1804 470 
## convN: -0.003951 m(z): 68 165 120 176 99 235 341 171 220 116 
## I: 200 
## convL: -0.224 n(z): 357 509 710 633 519 1103 3453 650 1913 403 
## convN: -0.006789 m(z): 67 166 119 177 99 235 339 171 222 116 
## I: 300 
## convL: -0.2263 n(z): 488 490 686 671 475 1078 3465 537 1996 364 
## convN: -0.003114 m(z): 67 168 119 177 99 235 340 171 219 116 
## I: 400 
## convL: -0.2382 n(z): 490 427 736 625 523 1122 3436 534 1945 412 
## convN: -0.003433 m(z): 68 167 122 175 99 234 339 171 220 116 
## I: 500 
## convL: -0.2215 n(z): 522 417 676 653 465 1150 3452 527 1999 389 
## convN: -0.00207 m(z): 67 168 120 176 98 237 337 171 221 116 
## I: 600 
## convL: -0.2197 n(z): 487 506 654 666 495 1168 3474 531 1895 374 
## convN: -0.005117 m(z): 68 167 121 176 99 234 341 170 218 117 
## I: 700 
## convL: -0.2165 n(z): 486 518 675 705 509 1135 3414 524 1930 354 
## convN: -0.003377 m(z): 67 167 120 177 98 236 341 171 218 116 
## I: 800 
## convL: -0.2134 n(z): 528 506 698 652 492 1090 3464 542 1948 330 
## convN: -0.00515 m(z): 68 169 118 177 104 232 336 170 221 116 
## Sample iterations: 200 
## I: 810 
## convL: -0.2304 n(z): 501 476 652 657 524 1123 3467 565 1908 377 
## convN: -0.004263 m(z): 68 167 119 177 99 236 338 171 220 116 
## I: 820 
## convL: -0.222 n(z): 513 500 759 662 469 1130 3464 570 1782 401 
## convN: -0.004598 m(z): 67 166 119 177 98 237 340 171 220 116 
## I: 830 
## convL: -0.2153 n(z): 514 516 699 711 461 1086 3425 532 1907 399 
## convN: -0.004271 m(z): 69 167 119 177 97 236 339 171 219 117 
## I: 840 
## convL: -0.216 n(z): 505 473 680 648 415 1222 3432 546 1881 448 
## convN: -0.002775 m(z): 69 168 119 177 97 236 339 170 218 118 
## I: 850 
## convL: -0.218 n(z): 500 442 723 626 448 1185 3554 513 1843 416 
## convN: -0.004411 m(z): 67 167 119 177 98 237 337 171 221 117 
## I: 860 
## convL: -0.2168 n(z): 524 490 667 618 452 1213 3472 536 1853 425 
## convN: -0.001134 m(z): 67 168 119 177 98 236 337 169 221 119 
## I: 870 
## convL: -0.2159 n(z): 574 466 680 592 430 1128 3458 557 1961 404 
## convN: -0.002191 m(z): 69 166 119 177 96 236 339 171 221 117 
## I: 880 
## convL: -0.2144 n(z): 500 505 678 596 485 1100 3429 549 1999 409 
## convN: -0.002394 m(z): 68 168 120 177 98 235 337 171 221 116 
## I: 890 
## convL: -0.2259 n(z): 530 511 708 653 458 1067 3415 510 2006 392 
## convN: -0.001423 m(z): 68 167 120 177 99 234 338 171 221 116 
## I: 900 
## convL: -0.2141 n(z): 506 515 685 589 447 1119 3513 531 1919 426 
## convN: -0.001895 m(z): 67 168 119 177 98 236 336 170 221 119 
## I: 910 
## convL: -0.2237 n(z): 538 509 672 667 443 1033 3468 548 1959 413 
## convN: -0.003129 m(z): 67 168 120 176 98 236 337 171 222 116 
## I: 920 
## convL: -0.2068 n(z): 562 507 681 672 411 1101 3471 536 1904 405 
## convN: -0.004273 m(z): 68 167 120 177 97 236 338 171 221 116 
## I: 930 
## convL: -0.2187 n(z): 517 501 651 643 444 1121 3447 510 2004 412 
## convN: -0.002495 m(z): 68 168 118 178 104 231 336 167 221 120 
## I: 940 
## convL: -0.2102 n(z): 529 480 652 556 481 1139 3475 562 1973 403 
## convN: -0.00344 m(z): 68 167 120 177 99 233 338 171 222 116 
## I: 950 
## convL: -0.2238 n(z): 522 479 652 600 478 1125 3441 531 1993 429 
## convN: -0.001377 m(z): 67 168 120 177 99 234 337 171 222 116 
## I: 960 
## convL: -0.2184 n(z): 547 482 642 597 474 1128 3463 528 1966 423 
## convN: -0.00399 m(z): 67 168 119 177 100 236 337 171 220 116 
## I: 970 
## convL: -0.2173 n(z): 540 474 676 580 459 1127 3424 522 1992 456 
## convN: -0.008498 m(z): 68 169 119 177 100 233 335 171 223 116 
## I: 980 
## convL: -0.2218 n(z): 563 446 663 561 440 1117 3453 541 2078 388 
## convN: -0.00307 m(z): 67 168 119 177 97 237 338 171 220 117 
## I: 990 
## convL: -0.2162 n(z): 521 461 658 592 476 1125 3519 554 1889 455 
## convN: -0.004939 m(z): 69 167 119 178 97 236 337 168 221 119 
## I: 1000 
## convL: -0.2198 n(z): 522 452 678 588 437 1182 3490 561 1916 424 
## convN: -0.001144 m(z): 68 168 120 177 99 234 337 171 221 116 
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
## [1] ggplot2_1.0.0       knitr_1.6           netresponse_1.17.11
## [4] reshape_0.8.5       mclust_4.4          minet_3.20.1       
## [7] igraph_0.7.1        Rgraphviz_2.8.1     graph_1.42.0       
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

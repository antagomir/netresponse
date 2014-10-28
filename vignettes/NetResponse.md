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
res <- generate.toydata()
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
```

![plot of chunk NetResponse2](figure/NetResponse2.png) 

### Network visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "network")
```

![plot of chunk NetResponse3](figure/NetResponse3.png) 

### Heatmap visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "heatmap")
```

![plot of chunk NetResponse4](figure/NetResponse4.png) 

### Boxplot visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "boxplot.data")
```

![plot of chunk NetResponse5](figure/NetResponse5.png) 

### Barplot visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "response.barplot")
```

```
## Error: $ operator is invalid for atomic vectors
```

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
subnet.id <- 'Subnet-2'
sample.probs <- response2sample(model, subnet.id)
```

```
## Error: object 'model' not found
```

```r
response.probs <- sample2response(model, subnet.id)
```

```
## Error: error in evaluating the argument 'model' in selecting a method for function 'getqofz': Error: object 'model' not found
```

Retrieve model parameters for a given subnetwork (Gaussian mixture
means, covariance diagonal, and component weights):


```r
get.model.parameters(model, subnet.id) 
```

```
## Error: object 'model' not found
```

## Extending the subnetworks

After identifying the locally connected subnetworks, it is possible to
search for features (genes) that are similar to a given subnetwork but
not directly interacting with it. To order the remaining features
in the input data based on similarity with the subnetwork, type


```r
g <- find.similar.features(model, subnet.id = "Subnet-1")
```

```
## Error: object 'model' not found
```

```r
subset(g, delta < 0)
```

```
## Error: object 'g' not found
```

This gives a data frame which indicates similarity level with the
subnetwork for each feature. The smaller, the more similar. Negative
values of delta indicate the presence of coordinated responses,
positive values of delta indicate independent responses. The data
frame is ordered such that the features are listed by decreasing
similarity.


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
## n(z): 1028 1019 992 1029 1019 1041 997 1069 1015 1041 
## m(z): 171 183 175 159 177 181 165 165 174 161 
## I: 100 
## convL: -0.2281 n(z): 480 360 805 568 3500 438 1310 541 516 1732 
## convN: -0.0018 m(z): 89 65 192 153 389 97 194 179 138 215 
## I: 200 
## convL: -0.2373 n(z): 500 350 787 567 3589 382 1033 577 488 1977 
## convN: -0.002454 m(z): 90 65 192 152 390 97 194 180 138 213 
## I: 300 
## convL: -0.2392 n(z): 500 385 835 523 3482 355 1084 619 535 1932 
## convN: -0.001107 m(z): 90 65 192 153 389 97 192 179 138 216 
## I: 400 
## convL: -0.2219 n(z): 455 366 851 506 3524 345 1062 599 508 2034 
## convN: -0.002304 m(z): 90 65 192 154 386 97 194 178 138 217 
## I: 500 
## convL: -0.2223 n(z): 462 372 812 546 3537 349 1080 612 514 1966 
## convN: -0.0005975 m(z): 90 65 192 154 389 97 194 178 138 214 
## I: 600 
## convL: -0.2159 n(z): 446 369 764 585 3524 363 1062 604 506 2027 
## convN: -0.001126 m(z): 88 65 192 153 388 97 194 180 138 216 
## I: 700 
## convL: -0.2218 n(z): 467 387 842 530 3484 380 1074 583 527 1976 
## convN: -0.001938 m(z): 88 65 192 152 390 97 195 181 138 213 
## I: 800 
## convL: -0.2268 n(z): 465 335 827 501 3498 366 1060 662 511 2025 
## convN: -0.002292 m(z): 88 65 192 154 389 97 194 179 138 215 
## Sample iterations: 200 
## I: 810 
## convL: -0.2329 n(z): 467 332 826 515 3526 319 1149 574 528 2014 
## convN: -0.001604 m(z): 88 65 192 153 391 97 194 180 138 213 
## I: 820 
## convL: -0.2303 n(z): 472 391 880 441 3552 338 1128 587 500 1961 
## convN: -0.0005098 m(z): 89 65 192 153 389 97 194 180 138 214 
## I: 830 
## convL: -0.2258 n(z): 491 363 808 499 3482 337 1098 601 480 2091 
## convN: -0.001747 m(z): 89 65 192 154 388 98 193 179 138 215 
## I: 840 
## convL: -0.2119 n(z): 443 348 857 481 3567 385 1045 617 526 1981 
## convN: -0.001167 m(z): 89 65 192 153 389 97 193 180 138 215 
## I: 850 
## convL: -0.2273 n(z): 474 354 832 428 3543 365 1180 651 523 1900 
## convN: -0.00132 m(z): 89 65 192 153 389 97 197 180 138 211 
## I: 860 
## convL: -0.2308 n(z): 431 353 867 472 3541 367 1133 582 518 1986 
## convN: -0.00235 m(z): 89 65 192 155 386 97 198 178 138 213 
## I: 870 
## convL: -0.2232 n(z): 486 354 851 510 3409 336 1120 617 516 2051 
## convN: -0.00115 m(z): 89 65 192 154 390 97 193 179 138 214 
## I: 880 
## convL: -0.2136 n(z): 441 362 847 480 3449 346 1054 664 505 2102 
## convN: -0.001212 m(z): 89 65 192 155 388 97 194 178 138 215 
## I: 890 
## convL: -0.2301 n(z): 470 352 825 470 3576 341 1095 610 534 1977 
## convN: -0.0004821 m(z): 89 65 192 153 385 97 201 180 138 211 
## I: 900 
## convL: -0.2287 n(z): 492 349 854 487 3576 378 1065 589 529 1931 
## convN: -0.005705 m(z): 88 65 191 155 388 97 197 178 138 214 
## I: 910 
## convL: -0.2278 n(z): 470 366 810 460 3575 365 1087 642 522 1953 
## convN: -0.003081 m(z): 89 65 191 154 389 97 193 179 138 216 
## I: 920 
## convL: -0.2175 n(z): 491 359 826 436 3573 380 1113 610 523 1939 
## convN: -0.002115 m(z): 89 65 192 152 389 97 195 181 138 213 
## I: 930 
## convL: -0.2259 n(z): 489 365 823 526 3528 318 1146 541 523 1991 
## convN: -0.002282 m(z): 89 65 192 153 391 97 193 180 138 213 
## I: 940 
## convL: -0.2277 n(z): 483 321 790 494 3618 326 1120 590 539 1969 
## convN: -0.003078 m(z): 89 65 191 153 389 97 193 181 138 215 
## I: 950 
## convL: -0.2146 n(z): 478 297 808 482 3614 358 1040 641 551 1981 
## convN: -0.0004351 m(z): 88 65 192 154 390 97 194 179 138 214 
## I: 960 
## convL: -0.2113 n(z): 478 364 801 447 3554 348 1063 631 528 2036 
## convN: -0.001617 m(z): 89 65 192 155 389 97 194 178 138 214 
## I: 970 
## convL: -0.2298 n(z): 466 369 785 465 3599 425 1073 638 525 1905 
## convN: -0.001457 m(z): 89 65 192 153 390 97 194 180 138 213 
## I: 980 
## convL: -0.2129 n(z): 451 422 766 458 3564 336 1102 678 520 1953 
## convN: -0.004158 m(z): 89 65 192 155 389 96 192 178 138 217 
## I: 990 
## convL: -0.2208 n(z): 418 352 785 468 3525 339 1167 632 535 2029 
## convN: -0.001885 m(z): 88 65 192 153 390 97 194 180 138 214 
## I: 1000 
## convL: -0.2287 n(z): 485 350 796 495 3520 309 1112 683 532 1968 
## convN: -0.001549 m(z): 89 65 191 154 387 97 196 179 138 215 
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
## [1] netresponse_1.17.1 reshape_0.8.5      mclust_4.4        
## [4] minet_3.20.1       igraph_0.7.1       Rgraphviz_2.8.1   
## [7] graph_1.42.0       knitr_1.6         
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

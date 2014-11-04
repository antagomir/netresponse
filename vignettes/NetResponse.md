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
## convL: -0.2356 n(z): 647 671 578 649 1218 555 1880 273 466 3313 
## convN: -0.001353 m(z): 119 180 175 171 238 122 202 68 101 335 
## I: 200 
## convL: -0.2211 n(z): 627 682 621 601 1145 602 1885 239 418 3430 
## convN: -0.001811 m(z): 120 182 177 170 239 122 204 69 95 333 
## I: 300 
## convL: -0.2397 n(z): 660 653 667 538 1149 555 1832 268 492 3436 
## convN: -0.003871 m(z): 120 175 176 171 235 122 204 67 98 343 
## I: 400 
## convL: -0.2276 n(z): 626 647 709 537 1159 558 1842 251 465 3456 
## convN: -0.002288 m(z): 118 174 176 172 235 121 206 67 97 345 
## I: 500 
## convL: -0.2177 n(z): 600 682 643 600 1106 558 1836 261 445 3519 
## convN: -0.003599 m(z): 120 175 177 170 237 121 203 67 97 344 
## I: 600 
## convL: -0.2207 n(z): 629 700 674 556 1068 553 1897 251 481 3441 
## convN: -0.004431 m(z): 119 174 177 170 232 121 207 68 98 345 
## I: 700 
## convL: -0.2031 n(z): 622 764 585 549 1016 620 1920 230 448 3496 
## convN: -0.00132 m(z): 120 174 177 170 236 121 204 67 96 346 
## I: 800 
## convL: -0.2198 n(z): 645 719 666 567 1035 600 1891 238 459 3430 
## convN: -0.001673 m(z): 120 174 176 171 236 121 204 67 96 346 
## Sample iterations: 200 
## I: 810 
## convL: -0.2215 n(z): 657 730 694 557 1047 579 1870 251 449 3416 
## convN: -0.004461 m(z): 117 174 177 171 230 121 208 68 100 345 
## I: 820 
## convL: -0.2197 n(z): 611 758 766 559 951 552 1947 265 452 3389 
## convN: -0.002284 m(z): 120 173 173 172 230 123 215 67 97 341 
## I: 830 
## convL: -0.2228 n(z): 613 698 708 568 1097 534 1975 256 475 3326 
## convN: -0.002003 m(z): 120 174 176 171 234 121 207 67 95 346 
## I: 840 
## convL: -0.2077 n(z): 654 733 703 517 995 545 2128 230 426 3319 
## convN: -0.002554 m(z): 120 174 177 170 232 121 210 67 95 345 
## I: 850 
## convL: -0.2116 n(z): 655 717 776 524 1003 572 1996 216 443 3348 
## convN: -0.007814 m(z): 118 175 176 172 232 120 210 67 98 343 
## I: 860 
## convL: -0.2203 n(z): 624 744 718 540 991 578 1933 237 471 3414 
## convN: -0.001764 m(z): 120 174 177 170 231 121 210 67 96 345 
## I: 870 
## convL: -0.2206 n(z): 610 746 724 517 1058 569 1961 247 462 3356 
## convN: -0.0022 m(z): 120 174 177 170 232 121 209 67 97 344 
## I: 880 
## convL: -0.2198 n(z): 584 774 743 482 1019 578 1949 217 447 3457 
## convN: -0.002698 m(z): 120 175 178 169 233 121 207 68 96 344 
## I: 890 
## convL: -0.2176 n(z): 626 722 707 522 1073 567 1933 241 419 3440 
## convN: -0.001619 m(z): 120 174 174 172 236 122 204 67 96 346 
## I: 900 
## convL: -0.2151 n(z): 617 747 689 575 1067 616 1796 225 443 3475 
## convN: -0.004143 m(z): 120 173 175 172 236 122 207 67 96 343 
## I: 910 
## convL: -0.2226 n(z): 643 796 701 537 1025 561 1876 269 455 3387 
## convN: -0.004811 m(z): 120 174 177 170 235 122 206 67 98 342 
## I: 920 
## convL: -0.2194 n(z): 641 751 695 537 1106 557 1858 246 486 3373 
## convN: -0.004781 m(z): 121 174 176 170 234 121 208 67 95 345 
## I: 930 
## convL: -0.2209 n(z): 618 762 635 522 1046 572 1951 219 449 3476 
## convN: -0.006343 m(z): 119 176 176 172 232 121 210 67 97 341 
## I: 940 
## convL: -0.2156 n(z): 616 722 669 551 981 602 1969 241 470 3429 
## convN: -0.004182 m(z): 120 174 177 170 234 122 209 67 96 342 
## I: 950 
## convL: -0.2297 n(z): 604 700 636 530 1068 565 1965 259 468 3455 
## convN: -0.003164 m(z): 120 176 177 170 236 121 205 68 97 341 
## I: 960 
## convL: -0.221 n(z): 596 766 643 530 1009 547 1949 255 449 3506 
## convN: -0.0025 m(z): 119 174 177 171 233 121 207 67 96 346 
## I: 970 
## convL: -0.2249 n(z): 605 854 629 561 985 567 1930 256 435 3428 
## convN: -0.003478 m(z): 120 175 175 172 232 121 208 67 97 344 
## I: 980 
## convL: -0.2179 n(z): 622 826 700 544 1023 567 1888 265 430 3385 
## convN: -0.004325 m(z): 120 178 177 170 233 122 206 67 97 341 
## I: 990 
## convL: -0.2106 n(z): 651 795 628 560 975 573 1880 274 459 3455 
## convN: -0.001699 m(z): 118 174 177 171 232 121 208 67 97 346 
## I: 1000 
## convL: -0.2173 n(z): 637 769 738 554 1044 535 1860 234 446 3433 
## convN: -0.001898 m(z): 120 175 177 170 234 121 208 67 96 343 
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
## [1] knitr_1.6           ggplot2_1.0.0       netresponse_1.17.12
## [4] reshape_0.8.5       mclust_4.4          minet_3.20.1       
## [7] Rgraphviz_2.8.1     graph_1.42.0       
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

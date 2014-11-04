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
## convL: -0.2553 n(z): 535 696 443 706 464 1843 3011 686 1102 764 
## convN: -0.00337 m(z): 67 175 121 177 99 237 343 170 203 119 
## I: 200 
## convL: -0.2528 n(z): 494 703 514 698 503 1870 2986 632 1061 789 
## convN: -0.004998 m(z): 67 176 120 177 101 238 343 171 201 117 
## I: 300 
## convL: -0.2286 n(z): 496 709 471 685 547 1891 3083 677 959 732 
## convN: -0.001944 m(z): 67 174 120 177 100 240 343 171 202 117 
## I: 400 
## convL: -0.2405 n(z): 487 689 447 669 604 1982 2930 652 955 835 
## convN: -0.003218 m(z): 67 176 120 178 107 238 339 169 200 117 
## I: 500 
## convL: -0.2336 n(z): 511 726 464 655 627 1989 3007 586 934 751 
## convN: -0.001557 m(z): 67 173 120 177 107 237 341 169 202 118 
## I: 600 
## convL: -0.2468 n(z): 506 659 522 700 625 1957 2939 603 932 807 
## convN: -0.004944 m(z): 67 174 123 175 105 237 342 168 201 119 
## I: 700 
## convL: -0.2385 n(z): 499 669 494 696 609 1893 3011 599 989 791 
## convN: -0.004978 m(z): 67 179 121 173 105 236 340 171 202 117 
## I: 800 
## convL: -0.2465 n(z): 545 681 486 694 584 1936 3026 618 947 733 
## convN: -0.003516 m(z): 67 176 123 173 105 237 340 170 202 118 
## Sample iterations: 200 
## I: 810 
## convL: -0.2415 n(z): 509 680 480 707 600 1962 2960 624 990 738 
## convN: -0.0067 m(z): 67 176 123 176 103 238 339 168 202 119 
## I: 820 
## convL: -0.2417 n(z): 477 698 478 699 602 1928 2944 654 1000 770 
## convN: -0.001836 m(z): 67 179 122 174 105 237 339 168 201 119 
## I: 830 
## convL: -0.2377 n(z): 466 654 516 670 585 1880 2991 670 1005 813 
## convN: -0.001716 m(z): 67 176 123 174 105 238 341 169 200 118 
## I: 840 
## convL: -0.2457 n(z): 478 727 485 639 616 1941 2981 704 925 754 
## convN: -0.004032 m(z): 67 179 122 174 105 237 341 168 199 119 
## I: 850 
## convL: -0.2375 n(z): 463 649 507 731 538 2055 3014 669 901 723 
## convN: -0.002892 m(z): 67 176 123 173 105 237 341 171 201 117 
## I: 860 
## convL: -0.2391 n(z): 449 657 517 699 579 1917 3020 667 911 834 
## convN: -0.00103 m(z): 67 175 124 173 105 237 342 172 200 116 
## I: 870 
## convL: -0.2314 n(z): 435 678 564 696 622 1897 2998 658 939 763 
## convN: -0.000995 m(z): 68 177 124 175 104 238 335 169 203 118 
## I: 880 
## convL: -0.2375 n(z): 465 704 481 697 584 1926 2967 642 973 811 
## convN: -0.004499 m(z): 67 178 122 175 104 237 340 169 201 118 
## I: 890 
## convL: -0.2448 n(z): 423 730 481 720 585 1988 3007 653 968 695 
## convN: -0.003068 m(z): 67 176 122 173 105 238 340 172 202 116 
## I: 900 
## convL: -0.2456 n(z): 459 698 513 679 570 1948 3017 653 982 731 
## convN: -0.004177 m(z): 67 177 123 173 105 237 341 172 200 116 
## I: 910 
## convL: -0.2308 n(z): 440 664 486 654 618 1980 2904 671 996 837 
## convN: -0.008117 m(z): 68 176 123 173 105 236 338 170 204 118 
## I: 920 
## convL: -0.2367 n(z): 460 724 472 654 650 1981 2930 664 895 820 
## convN: -0.001491 m(z): 67 176 123 173 105 237 339 172 203 116 
## I: 930 
## convL: -0.2297 n(z): 438 691 455 669 541 1972 3003 698 941 842 
## convN: -0.003354 m(z): 67 176 124 173 105 238 339 172 201 116 
## I: 940 
## convL: -0.2372 n(z): 456 724 472 672 596 2008 3000 659 934 729 
## convN: -0.003375 m(z): 67 177 122 173 105 237 339 171 203 117 
## I: 950 
## convL: -0.2431 n(z): 477 683 469 665 532 1999 3038 682 985 720 
## convN: -0.003865 m(z): 67 175 123 174 105 237 340 169 203 118 
## I: 960 
## convL: -0.2395 n(z): 462 669 457 690 605 1861 3027 683 973 823 
## convN: -0.001337 m(z): 68 177 121 174 105 238 339 171 202 116 
## I: 970 
## convL: -0.2434 n(z): 471 688 479 669 600 1917 3055 668 962 741 
## convN: -0.001059 m(z): 67 176 123 173 105 237 339 172 203 116 
## I: 980 
## convL: -0.2343 n(z): 447 690 452 714 588 1983 3028 677 940 731 
## convN: -0.002408 m(z): 67 177 123 174 105 238 338 169 202 118 
## I: 990 
## convL: -0.2495 n(z): 469 689 455 746 587 1899 3012 710 960 723 
## convN: -0.002227 m(z): 67 178 123 176 103 236 339 169 202 118 
## I: 1000 
## convL: -0.2328 n(z): 453 716 487 674 650 1890 3019 655 921 785 
## convN: -0.00398 m(z): 67 179 121 176 103 238 339 169 201 118 
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

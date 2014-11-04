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

![plot of chunk NetResponse3](figure/NetResponse3.png) 

### Heatmap visualization


```r
vis <- plot.responses(res, subnet.id, plot.mode = "heatmap")
```

```
## Error: 'at' and 'labels' lengths differ, 200 != 1
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
## convL: -0.2515 n(z): 591 559 487 749 610 1153 2956 633 1826 686 
## convN: -0.003911 m(z): 68 169 124 176 103 230 337 168 217 119 
## I: 200 
## convL: -0.2568 n(z): 588 656 485 735 528 1058 2924 632 1877 767 
## convN: -0.001785 m(z): 68 170 123 176 103 230 337 169 217 118 
## I: 300 
## convL: -0.2371 n(z): 579 603 488 751 519 1137 2901 656 1839 777 
## convN: -0.004979 m(z): 68 177 124 176 103 236 337 168 203 119 
## I: 400 
## convL: -0.2397 n(z): 588 600 525 767 526 1182 2883 638 1734 807 
## convN: -0.001519 m(z): 68 181 124 176 103 235 333 168 204 119 
## I: 500 
## convL: -0.2425 n(z): 566 580 493 769 535 1260 2941 609 1728 769 
## convN: -0.000779 m(z): 68 181 124 176 103 235 333 169 204 118 
## I: 600 
## convL: -0.2384 n(z): 588 682 501 772 516 1087 3053 624 1713 714 
## convN: -0.003301 m(z): 68 181 124 176 103 236 334 168 202 119 
## I: 700 
## convL: -0.2462 n(z): 619 678 491 713 530 1114 2985 653 1670 797 
## convN: -0.001854 m(z): 68 180 124 176 103 235 334 170 204 117 
## I: 800 
## convL: -0.2531 n(z): 564 655 508 736 552 1126 2868 609 1846 786 
## convN: -0.002105 m(z): 68 182 124 176 103 235 333 168 204 118 
## Sample iterations: 200 
## I: 810 
## convL: -0.2402 n(z): 576 662 507 699 535 1192 2857 618 1775 829 
## convN: -0.002229 m(z): 68 179 124 176 102 235 335 169 204 119 
## I: 820 
## convL: -0.2355 n(z): 566 639 549 718 541 1157 2873 614 1749 844 
## convN: -0.003324 m(z): 68 179 122 176 103 235 338 168 203 119 
## I: 830 
## convL: -0.2331 n(z): 572 671 480 732 532 1192 2930 654 1725 762 
## convN: -0.002438 m(z): 68 180 123 176 103 234 336 168 204 119 
## I: 840 
## convL: -0.2499 n(z): 596 695 504 744 504 1141 2925 593 1761 787 
## convN: -0.001927 m(z): 68 179 124 176 103 233 337 168 204 119 
## I: 850 
## convL: -0.2411 n(z): 578 703 493 749 496 1158 2964 590 1700 819 
## convN: -0.001629 m(z): 68 181 124 176 103 235 334 168 203 119 
## I: 860 
## convL: -0.2551 n(z): 530 700 527 773 559 1149 2927 597 1703 785 
## convN: -0.002122 m(z): 68 182 124 176 103 234 334 169 203 118 
## I: 870 
## convL: -0.245 n(z): 604 691 490 735 505 1217 2851 577 1751 829 
## convN: -0.002248 m(z): 68 181 122 176 103 236 336 167 203 119 
## I: 880 
## convL: -0.243 n(z): 618 636 521 716 556 1169 2903 620 1750 761 
## convN: -0.002882 m(z): 68 178 124 176 103 233 339 168 204 118 
## I: 890 
## convL: -0.2548 n(z): 589 648 459 769 562 1143 2909 628 1786 757 
## convN: -0.004656 m(z): 68 179 124 176 103 235 335 169 204 118 
## I: 900 
## convL: -0.2487 n(z): 625 616 495 700 554 1168 2914 612 1741 825 
## convN: -0.002511 m(z): 68 180 124 176 103 235 335 168 204 118 
## I: 910 
## convL: -0.2434 n(z): 570 674 457 717 528 1217 2877 645 1775 790 
## convN: -0.001541 m(z): 68 182 124 176 103 235 332 168 204 119 
## I: 920 
## convL: -0.2363 n(z): 545 683 517 750 552 1187 2918 622 1677 799 
## convN: -0.00367 m(z): 68 181 126 175 103 237 333 167 202 119 
## I: 930 
## convL: -0.2421 n(z): 588 660 485 739 562 1113 2962 617 1766 758 
## convN: -0.00118 m(z): 68 182 124 176 103 235 332 168 204 119 
## I: 940 
## convL: -0.2545 n(z): 604 666 492 721 509 1157 2888 638 1785 790 
## convN: -0.0008165 m(z): 68 179 124 176 103 235 335 168 204 119 
## I: 950 
## convL: -0.2318 n(z): 603 699 476 737 562 1050 2881 673 1773 796 
## convN: -0.001694 m(z): 68 180 125 176 103 235 334 168 204 118 
## I: 960 
## convL: -0.2533 n(z): 599 731 454 757 480 1084 2940 655 1800 750 
## convN: -0.002807 m(z): 68 177 124 176 103 233 340 168 203 119 
## I: 970 
## convL: -0.2472 n(z): 645 735 512 750 506 1154 2940 589 1679 740 
## convN: -0.002724 m(z): 68 182 124 176 103 237 332 169 202 118 
## I: 980 
## convL: -0.249 n(z): 597 729 483 731 534 1143 2896 606 1722 809 
## convN: -0.000472 m(z): 68 181 124 176 103 237 333 168 202 119 
## I: 990 
## convL: -0.2377 n(z): 571 708 503 721 512 1116 2979 647 1747 746 
## convN: -0.00238 m(z): 68 179 123 176 103 234 339 168 203 118 
## I: 1000 
## convL: -0.2384 n(z): 574 713 480 771 545 1088 2895 631 1729 824 
## convN: -0.00237 m(z): 68 179 123 176 103 232 341 168 203 118 
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

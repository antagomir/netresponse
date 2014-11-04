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
## convL: -0.2526 n(z): 560 711 446 569 570 1775 3225 616 1128 650 
## convN: -0.002933 m(z): 67 178 119 176 99 238 342 171 204 117 
## I: 200 
## convL: -0.2381 n(z): 514 729 460 543 597 1857 3189 624 1096 641 
## convN: -0.006424 m(z): 67 178 121 174 99 239 341 172 203 117 
## I: 300 
## convL: -0.2312 n(z): 468 704 474 605 625 1907 3306 608 989 564 
## convN: -0.002759 m(z): 67 177 119 176 100 239 342 172 202 117 
## I: 400 
## convL: -0.2413 n(z): 432 740 548 640 639 1887 3226 609 982 547 
## convN: -0.006176 m(z): 67 176 123 174 105 237 341 172 200 116 
## I: 500 
## convL: -0.2431 n(z): 438 711 559 596 575 1980 3223 544 1048 576 
## convN: -0.004918 m(z): 68 178 124 175 103 236 339 171 201 116 
## I: 600 
## convL: -0.2373 n(z): 439 720 566 605 608 2035 3168 543 986 580 
## convN: -0.002369 m(z): 67 178 124 174 103 237 339 172 201 116 
## I: 700 
## convL: -0.2291 n(z): 465 761 502 637 595 1989 3182 561 1022 536 
## convN: -0.004354 m(z): 68 178 124 173 104 236 339 172 201 116 
## I: 800 
## convL: -0.2337 n(z): 458 709 525 587 611 2034 3190 566 1001 569 
## convN: -0.0032 m(z): 67 178 124 174 103 237 337 173 202 116 
## Sample iterations: 200 
## I: 810 
## convL: -0.2333 n(z): 467 695 548 551 637 2073 3231 563 904 581 
## convN: -0.002335 m(z): 68 177 123 174 103 236 340 173 201 116 
## I: 820 
## convL: -0.249 n(z): 443 705 579 598 624 1958 3220 542 980 601 
## convN: -0.003527 m(z): 67 175 124 172 105 237 340 173 202 116 
## I: 830 
## convL: -0.2356 n(z): 459 747 564 552 670 1982 3157 535 1000 584 
## convN: -0.001876 m(z): 67 179 122 174 103 237 340 172 201 116 
## I: 840 
## convL: -0.2317 n(z): 445 721 600 577 614 1964 3235 527 979 588 
## convN: -0.00206 m(z): 68 179 122 174 103 236 340 173 200 116 
## I: 850 
## convL: -0.2261 n(z): 427 717 576 608 679 1929 3180 525 1006 603 
## convN: -0.004911 m(z): 67 179 122 174 103 237 338 172 202 117 
## I: 860 
## convL: -0.2351 n(z): 433 741 550 602 642 1975 3290 535 953 529 
## convN: -0.004531 m(z): 67 177 122 175 103 236 342 171 202 116 
## I: 870 
## convL: -0.2299 n(z): 432 705 541 633 708 2028 3138 486 947 632 
## convN: -0.004395 m(z): 67 178 124 174 103 238 341 172 198 116 
## I: 880 
## convL: -0.2345 n(z): 453 725 579 599 706 2078 3113 492 927 578 
## convN: -0.004046 m(z): 67 177 123 174 103 238 341 172 200 116 
## I: 890 
## convL: -0.2192 n(z): 441 693 547 639 682 2039 3177 462 959 611 
## convN: -0.002862 m(z): 67 178 123 174 104 238 338 171 201 117 
## I: 900 
## convL: -0.2271 n(z): 437 711 562 554 655 2024 3276 505 952 574 
## convN: -0.006498 m(z): 67 178 124 173 104 238 338 173 200 116 
## I: 910 
## convL: -0.229 n(z): 426 725 536 631 681 1999 3204 496 953 599 
## convN: -0.004949 m(z): 67 177 123 173 104 238 339 173 201 116 
## I: 920 
## convL: -0.2326 n(z): 414 741 569 600 617 2048 3171 518 975 597 
## convN: -0.004856 m(z): 67 178 123 173 105 237 339 170 202 117 
## I: 930 
## convL: -0.2277 n(z): 487 751 569 560 582 1980 3305 563 891 562 
## convN: -0.001738 m(z): 67 177 124 174 103 238 339 173 200 116 
## I: 940 
## convL: -0.2343 n(z): 440 739 527 552 633 1940 3217 529 1064 609 
## convN: -0.00558 m(z): 67 178 123 174 102 238 337 174 202 116 
## I: 950 
## convL: -0.2359 n(z): 439 783 590 561 651 2003 3101 561 957 604 
## convN: -0.004692 m(z): 67 179 122 175 103 238 340 171 199 117 
## I: 960 
## convL: -0.2239 n(z): 475 686 568 617 662 1986 3153 529 974 600 
## convN: -0.00291 m(z): 67 178 123 175 104 237 339 171 200 117 
## I: 970 
## convL: -0.238 n(z): 461 701 544 625 637 1929 3160 577 977 639 
## convN: -0.008339 m(z): 67 174 124 174 103 236 343 173 201 116 
## I: 980 
## convL: -0.2373 n(z): 435 759 568 591 662 1833 3205 600 973 624 
## convN: -0.001904 m(z): 67 178 124 174 103 236 340 172 201 116 
## I: 990 
## convL: -0.2392 n(z): 464 683 589 553 663 1915 3191 560 1016 616 
## convN: -0.001156 m(z): 67 178 123 174 103 238 339 173 200 116 
## I: 1000 
## convL: -0.2221 n(z): 437 672 601 582 690 1920 3276 543 945 584 
## convN: -0.00658 m(z): 67 178 123 173 103 240 339 173 199 116 
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

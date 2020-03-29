---
title: "Introduction to the netresponse R package"
author: "Leo Lahti et al."
date: "2020-03-29"
output:
  BiocStyle::html_document:
    toc: true
    fig_caption: yes    
  rmarkdown::md_document:
    toc: true
  rmarkdown::pdf_document:
    toc: true    
vignette: >
  %\VignetteIndexEntry{microbiome R package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




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
investigate the structure

NetResponse detects and characterizes
subnetworks that exhibit context-specific activation patterns across
versatile collections of functional measurements, such as gene
expression data. The implementations are partially based on the
agglomerative independent variable group analysis
([AIVGA](http://www.sciencedirect.com/science/article/pii/S0925231208000659))
and variational Dirichlet process Gaussian mixture models (Kurihara
et al. 2007). The
tools are particularly useful for global exploratory analysis of
genome-wide interaction networks and versatile collections of gene
expression data, and in general to discover subnetworks with
alternative states.


## Usage examples

Examples on running NetResponse algorithm and visualizing the
results. The algorithm combines network and functional information to
detect coherent subnetworks that reveal distinct activation modes
across conditions. 


```r
library(netresponse)

# Generate simulated data
res <- generate.toydata(Dim = 3, Nc = 3, Ns = 200, sd0 = 3, rgam.shape = 1, rgam.scale = 1)

D <- res$data
component.means <- res$means
component.sds   <- res$sds
sample2comp     <- res$sample2comp

# Use fully connected network
network <- matrix(rep(1, 9), nrow = 3) 

# Fit NetResponse model
# Various network formats are supported, see help(detect.responses) for
# details. With large data sets, consider the 'speedup' option.
set.seed(4243)
res <- detect.responses(D, network, mixture.method = "vdp", pca.basis = TRUE)

# List subnets (each is a list of nodes)
subnet.id <- names(get.subnets(res))[[1]]
```

### PCA visualization


```r
library(ggplot2)
vis <- plot_responses(res, subnet.id, plot_mode = "pca")
```


```r
# Modify the resulting ggplot2 object to enhance visualization
p <- vis$p # Pick the ggplot2 object from results
p <- p + geom_point(size = 3) # Modify point size
print(p) # Plot
```

![plot of chunk NetResponse2b](fig/NetResponse2b-1.png)


### Network visualization


```r
vis <- plot_responses(res, subnet.id, plot_mode = "network")
```

![plot of chunk NetResponse3](fig/NetResponse3-1.png)

### Heatmap visualization


```r
vis <- plot_responses(res, subnet.id, plot_mode = "heatmap")
```

![plot of chunk NetResponse4](fig/NetResponse4-1.png)

### Boxplot visualization


```r
vis <- plot_responses(res, subnet.id, plot_mode = "boxplot_data")
```

![plot of chunk NetResponse5](fig/NetResponse5-1.png)

See also mode = "response.barplot" 


### Color scale


```r
plot_scale(vis$breaks, vis$palette, two.sided = TRUE)
```

![plot of chunk NetResponse7](fig/NetResponse7-1.png)


### Cluster assignments

The sample-response assignments from the mixture model are soft
ie. defined as continuous probabilities. Retrieve the hard clustering
ie. list of samples for each response, response for each sample, based
the highest probability:


```r
subnet.id <- 'Subnet-1'

# Sample - response probabilities (soft cluster assignment)
response.probs <- sample2response(res, subnet.id)
tail(round(response.probs, 6))
```

```
##            Mode-1 Mode-2
## Sample-195      1      0
## Sample-196      1      0
## Sample-197      1      0
## Sample-198      1      0
## Sample-199      1      0
## Sample-200      1      0
```

```r
# Sample - response hard assignments
hard.clusters <- response2sample(res, subnet.id)
print(hard.clusters)
```

```
## $`Mode-1`
##   [1] "Sample-1"   "Sample-2"   "Sample-6"   "Sample-7"   "Sample-9"  
##   [6] "Sample-10"  "Sample-11"  "Sample-16"  "Sample-17"  "Sample-19" 
##  [11] "Sample-20"  "Sample-21"  "Sample-23"  "Sample-24"  "Sample-25" 
##  [16] "Sample-26"  "Sample-27"  "Sample-29"  "Sample-30"  "Sample-34" 
##  [21] "Sample-35"  "Sample-36"  "Sample-37"  "Sample-39"  "Sample-42" 
##  [26] "Sample-43"  "Sample-46"  "Sample-47"  "Sample-48"  "Sample-49" 
##  [31] "Sample-50"  "Sample-51"  "Sample-53"  "Sample-54"  "Sample-55" 
##  [36] "Sample-56"  "Sample-57"  "Sample-58"  "Sample-59"  "Sample-60" 
##  [41] "Sample-61"  "Sample-62"  "Sample-63"  "Sample-64"  "Sample-65" 
##  [46] "Sample-66"  "Sample-67"  "Sample-68"  "Sample-70"  "Sample-71" 
##  [51] "Sample-72"  "Sample-73"  "Sample-75"  "Sample-77"  "Sample-81" 
##  [56] "Sample-82"  "Sample-84"  "Sample-86"  "Sample-87"  "Sample-89" 
##  [61] "Sample-91"  "Sample-92"  "Sample-93"  "Sample-95"  "Sample-97" 
##  [66] "Sample-100" "Sample-101" "Sample-102" "Sample-103" "Sample-104"
##  [71] "Sample-105" "Sample-106" "Sample-107" "Sample-109" "Sample-110"
##  [76] "Sample-111" "Sample-112" "Sample-113" "Sample-116" "Sample-119"
##  [81] "Sample-122" "Sample-123" "Sample-125" "Sample-126" "Sample-127"
##  [86] "Sample-128" "Sample-131" "Sample-132" "Sample-133" "Sample-134"
##  [91] "Sample-135" "Sample-136" "Sample-138" "Sample-139" "Sample-140"
##  [96] "Sample-142" "Sample-145" "Sample-146" "Sample-147" "Sample-148"
## [101] "Sample-149" "Sample-152" "Sample-153" "Sample-154" "Sample-156"
## [106] "Sample-157" "Sample-158" "Sample-159" "Sample-160" "Sample-161"
## [111] "Sample-162" "Sample-164" "Sample-165" "Sample-166" "Sample-167"
## [116] "Sample-169" "Sample-170" "Sample-172" "Sample-173" "Sample-174"
## [121] "Sample-177" "Sample-178" "Sample-179" "Sample-180" "Sample-181"
## [126] "Sample-183" "Sample-184" "Sample-187" "Sample-188" "Sample-189"
## [131] "Sample-190" "Sample-192" "Sample-196" "Sample-197" "Sample-199"
## 
## $`Mode-2`
##  [1] "Sample-3"   "Sample-4"   "Sample-5"   "Sample-8"   "Sample-12" 
##  [6] "Sample-13"  "Sample-14"  "Sample-15"  "Sample-18"  "Sample-22" 
## [11] "Sample-28"  "Sample-31"  "Sample-32"  "Sample-33"  "Sample-38" 
## [16] "Sample-40"  "Sample-41"  "Sample-44"  "Sample-45"  "Sample-52" 
## [21] "Sample-69"  "Sample-74"  "Sample-76"  "Sample-78"  "Sample-79" 
## [26] "Sample-80"  "Sample-83"  "Sample-85"  "Sample-88"  "Sample-90" 
## [31] "Sample-94"  "Sample-96"  "Sample-98"  "Sample-99"  "Sample-108"
## [36] "Sample-114" "Sample-115" "Sample-117" "Sample-118" "Sample-120"
## [41] "Sample-121" "Sample-124" "Sample-129" "Sample-130" "Sample-137"
## [46] "Sample-141" "Sample-143" "Sample-144" "Sample-150" "Sample-151"
## [51] "Sample-155" "Sample-163" "Sample-168" "Sample-171" "Sample-175"
## [56] "Sample-176" "Sample-182" "Sample-185" "Sample-186" "Sample-191"
## [61] "Sample-193" "Sample-194" "Sample-195" "Sample-198" "Sample-200"
```

Retrieve model parameters for a given subnetwork (Gaussian mixture
means, covariance diagonal, and component weights; see
help(get.model.parameters) for details):


```r
params <- get.model.parameters(res, subnet.id) 
names(params)
```

```
## [1] "mu"          "sd"          "w"           "free.energy" "Nparams"    
## [6] "qofz"        "nodes"
```

## Nonparametric Gaussian mixture models

Nonparametric Gaussian mixtures with variational Dirichlet processes
based on implementations by Kurihara et
al. (2007)
and [Honkela et
al.](http://www.sciencedirect.com/science/article/pii/S0925231208000659).


```r
# Generate 2-dimensional simulated data with 3 clusters
res <- generate.toydata(Dim = 2, Nc = 3, Ns = 200, sd0 = 3, rgam.shape = 1, rgam.scale = 1)

D <- res$data
real.means <- res$means
real.sds   <- res$sds
real.sample2comp     <- res$sample2comp

# Infinite Gaussian mixture model with       
# Variational Dirichlet Process approximation       
mixt <- vdp.mixt( D )
            
# Centroids of the detected Gaussian components       
estimated.means <- mixt$posterior$centroids

# The colors denote the known clusters
# The blue ball denotes the original (known) cluster centroids and
# the triangle denotes the estimated centroids
plot(D, col = real.sample2comp, pch = 1)
points(real.means, col = "blue", pch = 16, cex = 2)
points(estimated.means, col = "blue", pch = 17, cex = 2)
```

![plot of chunk vdp](fig/vdp-1.png)

```r
# Hard mixture component assignment for each sample
estimated.sample2comp <- apply(mixt$posterior$qOFz, 1, which.max)

# Compare known and estimated mixture components
# (note that cluster indices may have switched due to unidentifiability)
# nearly all samples have one-to-one match between the real and estimated 
# clusters
head(table(estimated.sample2comp, real.sample2comp))
```

```
##                      real.sample2comp
## estimated.sample2comp  1  2  3
##                     1  0 71  9
##                     2 65  0  2
##                     3  0  0 53
```

### Citing NetResponse

Please cite [Lahti et al. (2010)](http://bioinformatics.oxfordjournals.org/content/26/21/2713) with the package. 


```r
citation("netresponse")
```

```
## 
## To cite netrespose (algorithm and package) in publications use:
## 
##   Leo Lahti et al. Global modeling of transcriptional responses in
##   interaction networks Bioinformatics 26(21):2713--20, 2010.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Article{,
##     title = {Global modeling of transcriptional responses in interaction networks},
##     author = {Leo Lahti and Juha E.A. Knuuttila and Samuel Kaski},
##     journal = {Bioinformatics},
##     year = {2010},
##     volume = {26},
##     issue = {21},
##     pages = {2713--20},
##   }
## 
## For ICMg functionality, please cite additionally the references listed
## in help(ICMg.combined.sampler). Thanks for Olli-Pekka Huovilainen and
## Antonio Gusmao for contributions to the R/C implementation of the
## netresponse algorithm and Juuso Parkkinen for ICMg.
```

### Version information

This document was written using:


```r
sessionInfo()
```

```
## R version 3.6.3 Patched (2020-03-11 r78037)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 19.10
## 
## Matrix products: default
## BLAS:   /home/lei/bin/R-patched/lib/libRblas.so
## LAPACK: /home/lei/bin/R-patched/lib/libRlapack.so
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
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] ggplot2_3.3.0       knitr_1.28          netresponse_1.47.2 
##  [4] reshape2_1.4.3      mclust_5.4.5        minet_3.44.1       
##  [7] Rgraphviz_2.30.0    graph_1.64.0        BiocGenerics_0.32.0
## [10] devtools_2.2.2      usethis_1.5.1      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.4         mvtnorm_1.1-0      lattice_0.20-40    prettyunits_1.1.1 
##  [5] ps_1.3.2           assertthat_0.2.1   rprojroot_1.3-2    digest_0.6.25     
##  [9] R6_2.4.1           plyr_1.8.6         backports_1.1.5    stats4_3.6.3      
## [13] evaluate_0.14      highr_0.8          pillar_1.4.3       rlang_0.4.5.9000  
## [17] rstudioapi_0.11    callr_3.4.2        Matrix_1.2-18      qvalue_2.18.0     
## [21] labeling_0.3       desc_1.2.0         splines_3.6.3      stringr_1.4.0     
## [25] igraph_1.2.5       munsell_0.5.0      xfun_0.12          compiler_3.6.3    
## [29] pkgconfig_2.0.3    pkgbuild_1.0.6     tidyselect_1.0.0   tibble_2.1.3      
## [33] fansi_0.4.1        crayon_1.3.4       dplyr_0.8.99.9002  withr_2.1.2       
## [37] MASS_7.3-51.5      gtable_0.3.0       lifecycle_0.2.0    magrittr_1.5      
## [41] scales_1.1.0       cli_2.0.2          stringi_1.4.6      farver_2.0.3      
## [45] fs_1.3.2           remotes_2.1.1      testthat_2.3.2     ellipsis_0.3.0    
## [49] vctrs_0.2.99.9010  RColorBrewer_1.1-2 tools_3.6.3        dmt_0.8.20        
## [53] glue_1.3.2         purrr_0.3.3        processx_3.4.2     pkgload_1.0.2     
## [57] colorspace_1.4-1   sessioninfo_1.1.1  memoise_1.1.0
```

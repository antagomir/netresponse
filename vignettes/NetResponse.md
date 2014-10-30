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
get.model.parameters(res, subnet.id) 
```

```
## $mu
##        [,1]    [,2]      [,3]
## [1,] -7.628  0.9387  0.016593
## [2,]  6.708  1.8635 -0.023908
## [3,]  1.774 -3.7510  0.008401
## 
## $sd
##        [,1]   [,2]   [,3]
## [1,] 0.6103 0.5060 0.4644
## [2,] 0.6628 0.5279 0.5864
## [3,] 0.5078 0.6605 0.7228
## 
## $w
## [1] 0.375 0.355 0.270
## 
## $free.energy
##       [,1]
## [1,] 821.1
## 
## $Nparams
## [1] 21
## 
## $qofz
##              [,1]       [,2]       [,3]
##   [1,]  1.000e+00 8.542e-102  2.290e-88
##   [2,] 8.755e-119  1.000e+00  1.590e-34
##   [3,] 1.613e-113  1.000e+00  1.880e-37
##   [4,]  1.376e-66  1.081e-37  1.000e+00
##   [5,]  1.000e+00  5.301e-94  3.508e-76
##   [6,] 4.123e-135  1.000e+00  1.197e-49
##   [7,]  1.000e+00 4.113e-105  3.659e-87
##   [8,] 2.252e-126  1.000e+00  1.428e-38
##   [9,]  1.000e+00 4.035e-103  2.706e-84
##  [10,]  1.000e+00 5.607e-105  4.330e-89
##  [11,]  4.672e-69  1.484e-41  1.000e+00
##  [12,] 4.900e-132  1.000e+00  3.740e-41
##  [13,]  3.143e-61  1.142e-25  1.000e+00
##  [14,]  1.000e+00 2.040e-101  1.835e-92
##  [15,]  1.000e+00  1.949e-98  3.939e-82
##  [16,]  2.112e-72  5.585e-35  1.000e+00
##  [17,]  1.000e+00  8.150e-94  1.136e-72
##  [18,]  1.000e+00  1.155e-96  5.053e-80
##  [19,]  1.000e+00  1.930e-85  2.988e-68
##  [20,] 2.332e-105  1.000e+00  9.994e-23
##  [21,] 1.459e-123  1.000e+00  4.540e-40
##  [22,]  4.510e-69  6.665e-41  1.000e+00
##  [23,] 4.381e-110  1.000e+00  2.455e-32
##  [24,]  1.414e-70  1.197e-34  1.000e+00
##  [25,]  1.000e+00 4.308e-100  2.028e-80
##  [26,]  1.000e+00 4.914e-109  7.270e-90
##  [27,] 3.486e-136  1.000e+00  6.064e-48
##  [28,]  1.000e+00 3.720e-105  3.587e-91
##  [29,] 1.326e-123  1.000e+00  2.171e-37
##  [30,] 9.110e-114  1.000e+00  3.303e-32
##  [31,]  1.000e+00  3.561e-99  3.477e-85
##  [32,]  6.421e-74  2.241e-44  1.000e+00
##  [33,]  1.000e+00  1.492e-98  4.712e-86
##  [34,]  1.000e+00 2.396e-104  6.150e-91
##  [35,] 8.704e-118  1.000e+00  2.487e-36
##  [36,]  1.000e+00 8.138e-103  4.939e-85
##  [37,]  1.000e+00 1.099e-107  2.385e-92
##  [38,] 1.263e-120  1.000e+00  2.039e-38
##  [39,]  7.372e-80  1.482e-48  1.000e+00
##  [40,] 5.884e-113  1.000e+00  2.035e-35
##  [41,]  1.000e+00 3.213e-116 1.178e-103
##  [42,]  1.200e-66  2.639e-34  1.000e+00
##  [43,]  1.000e+00  3.692e-91  3.786e-72
##  [44,] 1.903e-136  1.000e+00  1.734e-48
##  [45,] 2.124e-126  1.000e+00  3.390e-42
##  [46,]  7.246e-75  1.528e-41  1.000e+00
##  [47,] 4.745e-120  1.000e+00  3.739e-38
##  [48,]  1.000e+00 4.920e-108  3.203e-89
##  [49,]  1.000e+00 1.680e-109  2.322e-94
##  [50,]  6.848e-98  1.000e+00  1.339e-27
##  [51,]  2.476e-68  6.214e-25  1.000e+00
##  [52,]  1.000e+00  3.578e-94  2.851e-78
##  [53,] 8.987e-116  1.000e+00  7.906e-36
##  [54,]  9.468e-96  1.000e+00  1.066e-25
##  [55,]  5.164e-70  5.894e-35  1.000e+00
##  [56,] 9.687e-110  1.000e+00  2.837e-30
##  [57,] 7.382e-133  1.000e+00  8.634e-42
##  [58,]  1.178e-59  4.183e-37  1.000e+00
##  [59,] 3.034e-118  1.000e+00  5.761e-39
##  [60,]  8.058e-77  7.278e-35  1.000e+00
##  [61,]  1.000e+00 4.070e-110  5.416e-95
##  [62,]  1.000e+00  1.615e-97  4.301e-80
##  [63,]  1.000e+00 6.047e-119 1.704e-106
##  [64,] 9.858e-110  1.000e+00  1.183e-29
##  [65,] 5.484e-125  1.000e+00  1.643e-38
##  [66,]  4.540e-68  2.560e-43  1.000e+00
##  [67,] 6.985e-131  1.000e+00  6.018e-40
##  [68,] 1.695e-116  1.000e+00  5.052e-34
##  [69,]  9.073e-68  6.733e-44  1.000e+00
##  [70,] 1.299e-106  1.000e+00  8.081e-28
##  [71,]  1.000e+00 1.280e-102  7.525e-84
##  [72,] 2.458e-114  1.000e+00  7.711e-34
##  [73,] 4.588e-117  1.000e+00  7.527e-38
##  [74,]  7.830e-70  4.827e-37  1.000e+00
##  [75,]  1.547e-78  1.746e-37  1.000e+00
##  [76,]  1.000e+00  3.245e-91  5.417e-75
##  [77,] 1.793e-129  1.000e+00  2.223e-42
##  [78,]  1.000e+00  7.109e-90  3.777e-73
##  [79,] 1.553e-127  1.000e+00  7.759e-40
##  [80,] 1.420e-124  1.000e+00  5.413e-41
##  [81,]  1.000e+00 4.273e-105  1.682e-94
##  [82,]  8.720e-63  2.616e-29  1.000e+00
##  [83,]  1.000e+00 5.560e-104  1.283e-86
##  [84,]  1.002e-68  3.771e-30  1.000e+00
##  [85,] 3.037e-108  1.000e+00  4.158e-31
##  [86,]  3.350e-64  2.697e-21  1.000e+00
##  [87,]  1.000e+00 1.237e-104  2.923e-88
##  [88,]  1.000e+00 1.664e-111  1.886e-97
##  [89,]  1.000e+00  4.762e-86  7.461e-70
##  [90,]  3.937e-57  6.692e-29  1.000e+00
##  [91,]  2.906e-79  4.173e-49  1.000e+00
##  [92,]  1.000e+00  1.504e-89  2.143e-70
##  [93,] 4.962e-119  1.000e+00  6.242e-32
##  [94,]  1.000e+00 5.738e-102  8.629e-88
##  [95,]  9.938e-70  1.001e-38  1.000e+00
##  [96,]  1.425e-65  2.498e-40  1.000e+00
##  [97,]  2.810e-88  1.000e+00  1.849e-19
##  [98,] 1.398e-129  1.000e+00  6.334e-45
##  [99,]  1.000e+00 1.054e-111  2.939e-93
## [100,]  1.039e-80  4.682e-32  1.000e+00
## [101,] 1.518e-110  1.000e+00  5.425e-32
## [102,]  1.000e+00  4.858e-86  3.810e-67
## [103,] 1.061e-131  1.000e+00  1.461e-43
## [104,]  1.000e+00 4.710e-104  2.991e-81
## [105,] 2.632e-120  1.000e+00  3.404e-36
## [106,]  5.262e-68  4.759e-27  1.000e+00
## [107,]  1.000e+00 2.368e-113  2.163e-97
## [108,]  1.802e-74  2.227e-37  1.000e+00
## [109,]  1.000e+00  6.161e-98  8.839e-83
## [110,] 1.076e-129  1.000e+00  3.972e-44
## [111,]  7.625e-69  1.292e-38  1.000e+00
## [112,]  9.660e-72  3.840e-37  1.000e+00
## [113,] 9.398e-119  1.000e+00  3.092e-34
## [114,] 9.389e-130  1.000e+00  3.624e-41
## [115,]  6.962e-64  2.106e-40  1.000e+00
## [116,]  1.000e+00  4.521e-90  3.414e-73
## [117,] 1.150e-113  1.000e+00  2.590e-36
## [118,]  1.000e+00  4.717e-96  1.191e-77
## [119,]  1.000e+00 1.171e-104  3.213e-85
## [120,] 3.549e-109  1.000e+00  4.481e-29
## [121,] 4.462e-122  1.000e+00  3.681e-39
## [122,]  4.335e-68  3.591e-25  1.000e+00
## [123,]  4.408e-61  2.090e-39  1.000e+00
## [124,]  1.000e+00  1.880e-98  3.613e-80
## [125,]  3.847e-60  3.344e-42  1.000e+00
## [126,]  1.000e+00 4.847e-106  8.525e-85
## [127,] 3.022e-124  1.000e+00  1.321e-38
## [128,]  1.000e+00  3.094e-95  6.028e-73
## [129,] 9.207e-100  1.000e+00  1.855e-25
## [130,]  1.000e+00 8.195e-105  5.043e-84
## [131,]  4.138e-65  2.135e-36  1.000e+00
## [132,]  1.000e+00  1.590e-83  2.675e-66
## [133,] 2.219e-121  1.000e+00  8.137e-32
## [134,] 1.608e-108  1.000e+00  6.363e-32
## [135,] 3.019e-121  1.000e+00  1.393e-37
## [136,]  6.863e-60  1.000e-30  1.000e+00
## [137,]  8.174e-76  3.137e-34  1.000e+00
## [138,]  1.000e+00 4.008e-100  1.713e-78
## [139,]  1.000e+00 2.901e-104  1.576e-88
## [140,]  1.000e+00  4.052e-91  9.870e-73
## [141,]  1.910e-78  3.207e-33  1.000e+00
## [142,]  1.000e+00 2.641e-103  1.016e-92
## [143,]  6.118e-96  1.000e+00  6.518e-20
## [144,]  1.000e+00 2.416e-102  3.378e-86
## [145,]  1.000e+00  5.316e-82  2.409e-64
## [146,]  1.000e+00  7.061e-92  7.076e-75
## [147,]  6.465e-76  2.531e-42  1.000e+00
## [148,]  1.000e+00 5.965e-107  3.538e-94
## [149,]  3.155e-60  3.967e-27  1.000e+00
## [150,]  3.874e-61  6.089e-40  1.000e+00
## [151,]  1.000e+00 4.885e-110 6.104e-101
## [152,]  1.000e+00 3.785e-105  1.128e-84
## [153,]  1.000e+00 9.254e-100  5.166e-84
## [154,]  3.822e-58  1.068e-38  1.000e+00
## [155,]  1.000e+00  1.169e-93  8.154e-78
## [156,]  1.000e+00  1.982e-82  2.070e-60
## [157,]  1.000e+00  2.405e-81  2.033e-58
## [158,]  8.221e-80  5.169e-42  1.000e+00
## [159,] 1.644e-117  1.000e+00  2.374e-34
## [160,]  1.561e-76  1.900e-40  1.000e+00
## [161,] 2.428e-116  1.000e+00  1.984e-31
## [162,]  1.000e+00  1.498e-92  3.651e-70
## [163,]  1.000e+00 2.985e-101  1.086e-85
## [164,]  5.065e-62  1.537e-33  1.000e+00
## [165,]  8.426e-80  4.382e-39  1.000e+00
## [166,] 1.013e-133  1.000e+00  4.381e-48
## [167,] 5.116e-126  1.000e+00  2.403e-41
## [168,] 4.071e-118  1.000e+00  4.613e-39
## [169,]  1.000e+00 3.230e-104  1.102e-87
## [170,]  1.000e+00 2.882e-107  1.271e-86
## [171,]  1.000e+00  6.978e-97  3.061e-80
## [172,] 2.300e-123  1.000e+00  1.342e-36
## [173,]  9.992e-70  5.589e-38  1.000e+00
## [174,]  1.701e-69  1.131e-26  1.000e+00
## [175,] 3.369e-120  1.000e+00  9.445e-38
## [176,] 1.312e-110  1.000e+00  3.768e-31
## [177,] 3.402e-101  1.000e+00  3.441e-26
## [178,]  6.037e-73  1.313e-37  1.000e+00
## [179,]  1.000e+00 4.026e-104  6.673e-91
## [180,]  1.334e-99  1.000e+00  1.313e-19
## [181,] 2.705e-100  1.000e+00  5.442e-29
## [182,] 1.074e-136  1.000e+00  2.278e-47
## [183,]  1.000e+00  9.951e-95  2.179e-70
## [184,]  1.000e+00 2.312e-115  6.698e-99
## [185,]  4.654e-62  4.088e-40  1.000e+00
## [186,] 5.340e-104  1.000e+00  4.704e-29
## [187,] 2.927e-126  1.000e+00  4.399e-39
## [188,]  3.860e-79  1.078e-33  1.000e+00
## [189,] 4.907e-121  1.000e+00  1.792e-35
## [190,]  1.000e+00 8.998e-101  5.879e-88
## [191,]  1.000e+00 3.076e-103  5.046e-83
## [192,]  1.000e+00  3.372e-99  2.788e-83
## [193,] 5.792e-130  1.000e+00  6.372e-44
## [194,]  1.000e+00 3.054e-108  1.304e-93
## [195,] 4.889e-124  1.000e+00  3.151e-39
## [196,]  3.893e-64  6.012e-52  1.000e+00
## [197,]  2.858e-66  8.378e-24  1.000e+00
## [198,]  2.447e-68  6.598e-36  1.000e+00
## [199,] 1.346e-127  1.000e+00  1.249e-40
## [200,] 1.376e-104  1.000e+00  2.528e-29
## 
## $nodes
## [1] "Feature-1" "Feature-2" "Feature-3"
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
## n(z): 1028 1019 992 1029 1019 1041 997 1069 1015 1041 
## m(z): 171 183 175 159 177 181 165 165 174 161 
## I: 100 
## convL: -0.246 n(z): 467 325 802 1388 3607 600 1011 728 478 844 
## convN: -0.005091 m(z): 110 69 169 162 362 116 205 184 95 239 
## I: 200 
## convL: -0.2346 n(z): 399 382 786 1378 3562 588 1099 746 488 822 
## convN: -0.0005401 m(z): 118 69 167 164 359 115 209 172 101 237 
## I: 300 
## convL: -0.2218 n(z): 408 385 826 1446 3576 630 1032 692 480 775 
## convN: -0.002696 m(z): 116 69 166 163 357 115 219 173 99 234 
## I: 400 
## convL: -0.2288 n(z): 431 389 834 1483 3538 593 1074 659 476 773 
## convN: -0.001766 m(z): 116 68 168 164 356 115 219 173 97 235 
## I: 500 
## convL: -0.2219 n(z): 418 399 824 1381 3554 642 1093 674 499 766 
## convN: -0.0009223 m(z): 116 67 168 164 357 115 218 173 97 236 
## I: 600 
## convL: -0.2377 n(z): 435 387 812 1281 3614 607 1170 660 499 785 
## convN: -0.001521 m(z): 116 67 166 163 360 116 219 173 98 233 
## I: 700 
## convL: -0.2225 n(z): 439 393 790 1361 3565 621 1194 599 502 786 
## convN: -0.003104 m(z): 116 68 166 164 359 115 218 173 98 234 
## I: 800 
## convL: -0.2216 n(z): 435 379 725 1437 3544 618 1157 676 466 813 
## convN: -0.004275 m(z): 116 67 168 164 358 115 219 173 97 234 
## Sample iterations: 200 
## I: 810 
## convL: -0.224 n(z): 453 377 783 1356 3469 648 1189 648 483 844 
## convN: -0.001356 m(z): 116 67 168 164 357 115 219 173 97 235 
## I: 820 
## convL: -0.2322 n(z): 435 334 800 1416 3511 642 1172 623 531 786 
## convN: -0.002259 m(z): 116 67 166 164 359 115 219 173 97 235 
## I: 830 
## convL: -0.2241 n(z): 437 352 821 1409 3598 617 1100 586 512 818 
## convN: -0.004005 m(z): 116 67 168 164 357 114 218 174 97 236 
## I: 840 
## convL: -0.2218 n(z): 399 330 784 1511 3618 617 1047 582 504 858 
## convN: -0.001349 m(z): 116 69 168 165 357 115 218 173 95 235 
## I: 850 
## convL: -0.2279 n(z): 374 343 786 1514 3568 648 1019 665 488 845 
## convN: -0.001582 m(z): 115 67 167 165 359 115 218 173 97 235 
## I: 860 
## convL: -0.23 n(z): 416 325 812 1486 3561 641 1080 545 481 903 
## convN: -0.001076 m(z): 116 67 167 165 358 115 218 173 96 236 
## I: 870 
## convL: -0.2345 n(z): 443 355 807 1436 3620 588 1033 641 454 873 
## convN: -0.001629 m(z): 116 69 167 165 358 116 218 172 95 235 
## I: 880 
## convL: -0.2186 n(z): 395 333 788 1488 3471 589 1086 690 507 903 
## convN: -0.001235 m(z): 115 68 167 165 358 115 218 173 98 234 
## I: 890 
## convL: -0.2306 n(z): 398 332 782 1411 3585 616 1073 688 455 910 
## convN: -0.001856 m(z): 115 67 168 165 357 115 218 173 97 236 
## I: 900 
## convL: -0.221 n(z): 392 336 747 1432 3604 619 1010 772 489 849 
## convN: -0.00526 m(z): 115 67 168 165 357 116 219 172 97 235 
## I: 910 
## convL: -0.2269 n(z): 440 343 742 1475 3531 582 1053 752 531 801 
## convN: -0.001782 m(z): 115 67 167 165 357 115 219 173 97 236 
## I: 920 
## convL: -0.2356 n(z): 433 369 830 1352 3536 574 1121 683 499 853 
## convN: -0.0009268 m(z): 116 69 168 165 357 115 218 173 95 235 
## I: 930 
## convL: -0.2155 n(z): 397 348 831 1439 3570 605 1070 697 475 818 
## convN: -0.002264 m(z): 115 67 166 165 360 115 219 173 97 234 
## I: 940 
## convL: -0.221 n(z): 401 364 810 1412 3458 581 1166 682 496 880 
## convN: -0.003276 m(z): 115 67 167 164 357 116 219 173 97 236 
## I: 950 
## convL: -0.2318 n(z): 405 351 810 1371 3478 622 1108 713 486 906 
## convN: -0.001478 m(z): 116 67 167 165 357 115 219 173 96 236 
## I: 960 
## convL: -0.2256 n(z): 394 353 833 1329 3546 600 1158 655 510 872 
## convN: -0.002093 m(z): 116 67 168 164 358 115 218 173 97 235 
## I: 970 
## convL: -0.2206 n(z): 443 326 810 1372 3557 626 1089 673 484 870 
## convN: -0.001074 m(z): 115 67 168 165 355 115 220 173 97 236 
## I: 980 
## convL: -0.2226 n(z): 428 345 834 1366 3453 642 1229 668 443 842 
## convN: -0.002021 m(z): 115 68 167 165 355 115 221 173 99 233 
## I: 990 
## convL: -0.229 n(z): 389 381 767 1379 3582 634 1030 652 542 894 
## convN: -0.002952 m(z): 115 67 166 165 357 116 220 172 98 235 
## I: 1000 
## convL: -0.2198 n(z): 435 351 845 1426 3546 612 1048 573 530 884 
## convN: -0.004654 m(z): 116 67 168 164 354 116 222 173 96 235 
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

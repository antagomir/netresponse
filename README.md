netresponse - Functional network analysis
===========

NetResponse provides a global map of network activation patterns. The
implementation is based on probabilistic models and variational
learning. Currently available tools:

 * netresponse
   ([1](http://bioinformatics.oxfordjournals.org/content/26/21/2713))
   to detect and characterize context-specific activation patterns and
   to construct global functional maps of large interaction networks
   by combining functional information with interaction networks

 * Interaction Component Model
   [ICMg](http://www.biomedcentral.com/1752-0509/4/4): to discover
   functional network modules, or communities, taking into account the
   uncertainty in network structure. Module discovery can be
   supervised by functional information of the network.

Applicability of these models has been demonstrated by case studies in
computational biology ([Lahti et
al. 2010](http://bioinformatics.oxfordjournals.org/content/26/21/2713),
[Parkkinen and Kaski
2010](http://www.biomedcentral.com/1752-0509/4/4), where the
algorithms have been used to investigate the structure and
context-specific transcriptional activity of genome-scale interaction
networks in human body.

The techniques are implemented in R. For installation instructions,
see the [project page at
BioConductor](http://www.bioconductor.org/help/bioc-views/devel/bioc/html/netresponse.html). For
further documentation, see the [package
vignette](http://www.bioconductor.org/packages/2.7/bioc/vignettes/netresponse/inst/doc/netresponse.pdf). In
addition, [Matlab
implementation](http://www.cis.hut.fi/projects/mi/software/NetResponse)
is available, but not supported.



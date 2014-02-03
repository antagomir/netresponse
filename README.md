netresponse 
===========

Tools for functional network analysis

netresponse provides a global map of network activation patterns. The
implementation is based on probabilistic models and variational
learning. Currently available tools include:

 * [netresponse] [Lahti et al. Bioinformatics
   2010](http://bioinformatics.oxfordjournals.org/content/26/21/2713))
   to detect and characterize context-specific activation patterns and
   to construct global functional maps of large interaction networks
   by combining functional information with interaction networks

 * Interaction Component Model ICMg [Parkkinen and Kaski. BMC Systems
   Biology 2010](http://www.biomedcentral.com/1752-0509/4/4): to
   discover functional network modules, or communities, taking into
   account the uncertainty in network structure. Module discovery can
   be supervised by functional information of the network.

Applicability of these models has been demonstrated by case studies in
computational biology (see the links above), where the algorithms have
been used to investigate the structure and context-specific
transcriptional activity of genome-scale interaction networks in human
body.

The techniques are implemented in R. For installation instructions,
see the [project page at
BioConductor](http://www.bioconductor.org/help/bioc-views/devel/bioc/html/netresponse.html). For
further documentation, see the [package
vignette](https://github.com/antagomir/netresponse/blob/master/vignettes/netresponse.pdf?raw=true). In
addition, [Matlab
implementation](http://www.cis.hut.fi/projects/mi/software/NetResponse)
is available, but not supported.


### Authors

Maintainer: [Leo Lahti](http://antagomir.github.io/info/contact)

Contributors: Olli-Pekka Huovilainen, Antonio Gusmao, Juuso Parkkinen



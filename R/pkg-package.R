#' @title Dna damage data set (PPI and expression) 
#' @description A combined yeast data set with protein-protein interactions and gene expression (dna damage). Gene expression profiles are transformed into links by computing a Pearson correlation for all pairs of genes and treating all correlations above 0.85 as additional links. Number of genes: 1823, number of interactions: 12382, number of gene expression observations: 52, number of total links with PPI and expression links: 15547.
#' @name dna
#' @docType data
#' @usage data(dna)
#' @format List of following objects: \describe{ \item{ppi}{PPI data matrix}
#' \item{exp}{gene expression profiles data matrix} \item{gids}{Vector of gene
#' ids corresponding to indices used in data matrices} \item{obs}{Gene
#' expression observation details} \item{combined.links}{pooled matrix of PPI
#' and expression links} }
#' @references Ulitsky, I. and Shamir, R. \emph{Identification of functional
#' modules using network topology and high-throughput data.} BMC Systems
#' Biology 2007, 1:8.
#' 
#' Nariai, N., Kolaczyk, E. D. and Kasif, S. \emph{Probabilistic Protein
#' Function Predition from Heterogenous Genome-Wide Data}. PLoS ONE 2007,
#' 2(3):e337.
#' 
#' Gasch, A., Huang, M., Metzner, S., Botstein, D. and Elledge, S.
#' \emph{Genomic expression responses to DNA-damaging agents and the regulatory
#' role of the yeast ATR homolog Mex1p.} Molecular Biology of the Cell 2001,
#' 12:2987-3003.
#' @source PPI data pooled from yeast data sets of [1] and [2]. Dna damage
#' expression set of [3].
#' @keywords datasets
#' @examples data(dna)
NULL


#' @title Class "NetResponseModel"
#' @description A NetResponse model.
#' @name NetResponseModel-class
#' @aliases NetResponseModel-class [[,NetResponseModel-method show,NetResponseModel-method 
#' @docType class
#' @section Objects from the Class: Returned by \code{\link{detect.responses}}
#' function.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @keywords classes
#' @examples showClass("NetResponseModel")
NULL


#' @title NetResponse: Global modeling of transcriptional responses in interaction networks
#' @description Global modeling of transcriptional responses in interaction networks.
#' \tabular{ll}{ Package: \tab netresponse\cr Type: \tab Package\cr Version:
#' \tab See sessionInfo() or DESCRIPTION file\cr Date: \tab 2011-02-03\cr
#' License: \tab GNU GPL >=2\cr LazyLoad: \tab yes\cr }
#' @name netresponse-package
#' @aliases netresponse-package netresponse
#' @docType package
#' @useDynLib netresponse
#' @author Leo Lahti, Olli-Pekka Huovilainen, Antonio Gusmao and Juuso
#' Parkkinen. Maintainer: Leo Lahti \email{leo.lahti@@iki.fi}
#' @references Leo Lahti et al.: Global modeling of transcriptional responses
#' in interaction networks. Bioinformatics (2010). See citation("netresponse")
#' for details.
#' @keywords package
#' @examples
#' # Define parameters for toy data
#' Ns <- 200  # number of samples (conditions)
#' Nf <- 10   # number of features (nodes)
#' feature.names <- paste("feat", seq(Nf), sep="")
#' sample.names  <- paste("sample", seq(Ns), sep="") 
#' # random seed
#' set.seed( 123 )
#' # Random network
#' netw <- pmax(array(sign(rnorm(Nf^2)), dim = c(Nf, Nf)), 0)
#' # in pathway analysis nodes correspond to genes
#' rownames(netw) <- colnames(netw) <- feature.names
#' # Random responses of the nodes across conditions 
#' D <- array(rnorm(Ns*Nf), dim = c(Ns,Nf), dimnames = list(sample.names, feature.names))
#' D[1:100, 4:6]  <- t(sapply(1:(Ns/2),function(x){rnorm(3, mean = 1:3)}))
#' D[101:Ns, 4:6] <- t(sapply(1:(Ns/2),function(x){rnorm(3, mean = 7:9)}))
#' # Calculate the model
#' model <- detect.responses(D, netw)
#' # Subnets (each is a list of nodes)
#' get.subnets( model )
#' 
#' # Retrieve model for one subnetwork
#' # means, standard devations and weights for the components
#' inds <- which(sapply(model@@last.grouping, length) > 2)
#' subnet.id <- names(model@@subnets)[[1]]
#' m <- get.model.parameters(model, subnet.id) 
#' print(m)
NULL

#' @title Osmoshock data set (PPI and expression)
#' @description A combined yeast data set with protein-protein interactions and gene expression (osmotick shock response). Gene expression profiles are transformed into links by computing a Pearson correlation for all pairs of genes and treating all correlations above 0.85 as additional links. Number of genes: 1711, number of interactions: 10250, number of gene expression observations: 133, number of total links with PPI and expression links: 14256.
#' @name osmo
#' @docType data
#' @usage data(osmo)
#' @format List of following objects: \describe{ \item{ppi}{PPI data matrix}
#' \item{exp}{gene expression profiles data matrix} \item{gids}{Vector of gene
#' ids corresponding to indices used in data matrices} \item{obs}{Gene
#' expression observation details} \item{combined.links}{pooled matrix of PPI
#' and expression links} }
#' @references Ulitsky, I. and Shamir, R. \emph{Identification of functional
#' modules using network topology and high-throughput data.} BMC Systems
#' Biology 2007, 1:8.
#' 
#' Nariai, N., Kolaczyk, E. D. and Kasif, S. \emph{Probabilistic Protein
#' Function Predition from Heterogenous Genome-Wide Data}. PLoS ONE 2007,
#' 2(3):e337.
#' 
#' O'Rourke, S. and Herskowitz, I. \emph{Unique and redundant roles for Hog
#' MAPK pathway components as revealed by whole-genome expression analysis.}
#' Molecular Biology of the Cell 2004, 15:532-42.
#' @source PPI data pooled from yeast data sets of [1] and [2]. Dna damage
#' expression set of [3].
#' @keywords datasets
#' @examples data(osmo)
NULL



#' @title toydata
#' @description Toy data for NetResponse examples.
#' @name toydata
#' @docType data
#' @usage data(toydata)
#' @format
#' Toy data: a list with three elements:
#' 
#' emat: Data matrix (samples x features). This contains the same features that
#' are provided in the network (toydata$netw). The matrix characterizes
#' measurements of network states across different conditions.
#' 
#' netw: Binary matrix that describes pairwise interactions between features.
#' This defines an undirected network over the features. A link between two
#' nodes is denoted by 1.
#' 
#' model: A pre-calculated model. Object of NetResponseModel class, resulting
#' from applying the netresponse algorithm on the toydata with model <-
#' detect.responses(D, netw).
#' @references Leo Lahti et al.: Global modeling of transcriptional responses
#' in interaction networks. Bioinformatics (2010).
#' @keywords misc
#' @examples
#'   data(toydata)
#'   D    <- toydata$emat   # Response matrix (samples x features)
#'   netw <- toydata$netw   # Network between the features
#'   model <- toydata$model # Pre-calculated NetResponseModel obtained with
#'                          # model <- detect.responses(D, netw)
NULL




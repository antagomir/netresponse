\name{ICMg.links.sampler}
\alias{ICMg.links.sampler}
\title{ ICMg: inferring gene modules from interaction data. }
\description{
  ICMg.links.sampler computes samples from the posterior of the assignments of datapoints (interactions) to latent components. From these we can then obtain component membership distributions and clusterings for genes.
}
\usage{
ICMg.links.sampler(L, C, alpha=10, beta=0.01,  B.num=8, B.size=100, S.num=20, S.size=10, C.boost=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{L}{ N x 2 matrix of link endpoints (N = number of links). }
  \item{C}{ Number of components. }
  \item{alpha}{ Hyperparameter describing the global distribution over components, larger alpha gives a more uniform distribution. }
  \item{beta}{ Hyperparameter describing the component-wise distributions over nodes, larger beta gives a more uniform distribution. }
  \item{B.num}{ Number of burnin rounds.* }
  \item{B.size}{ Size of one burnin round.* }
  \item{S.num}{ Number of sample rounds.* }
  \item{S.size}{ Size of one sample round.* }
  \item{C.boost}{ Set to 1 to use faster iteration with C, set to 0 to use slower R functions. }
}
\details{
  * One run consists of two parts, during burnin the sampler is expected to mix, after which the samples are taken. Information about convergence (convN and convL are estimates of convergence for link and node sampling, respectively) and component sizes are printed after each burnin/sample round. For example: B.num=8, B.size=100, S.num=20, S.size=10, runs 800 burnin iterations in 8 rounds and then takes 20 samples with an interval of 10 iterations.
}
\value{
  Returns samples as a list:
  \item{z}{ S.num x N matrix of samples of component assignments for links. }
  \item{conv}{ Vector of length (B.num + S.num) with convergence estimator values for link sampling. } 
  \item{counts}{ (B.num + S.num) x C matrix of link component sizes. } 
  additionally all parameters of the run are included in the list.
}
\author{ Juuso Parkkinen }

\references{
Parkkinen, J. and Kaski, S. Searching for functional gene modules with interaction component models.
BMC Systems Biology 4 (2010), 4.
}

\seealso{ \code{\link{ICMg.combined.sampler}}}
\examples{
	## Load data and set parameters
	data(osmo)
	C.boost = 1 ## Use faster C funtions
	alpha = 10
	beta = 0.01
	B.num = 10
	B.size = 10
	S.num = 10  
	S.size = 10
 	C = 24

	## Run sampling
	res = ICMg.links.sampler(osmo$ppi, C, alpha, beta, B.num, B.size, S.num, S.size, C.boost) 
	## Compute component membership probabilities for nodes
  	res$comp.memb <- ICMg.get.comp.memberships(osmo$ppi, res)
  	## Compute (hard) clustering for nodes
  	res$clustering <- apply(res$comp.memb, 2, which.max)
}
\keyword{ methods }

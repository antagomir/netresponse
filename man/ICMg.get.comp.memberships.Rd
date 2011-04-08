\name{ICMg.get.comp.memberships}
\alias{ICMg.get.comp.memberships}
\title{ ICMg.get.comp.memberships }

\description{Function for computing the component memberships for each data point from the MCMC samples.}
\usage{ ICMg.get.comp.memberships(links, samples) }

\author{ Juuso Parkkinen }

\references{
Parkkinen, J. and Kaski, S. Searching for functional gene modules with interaction component models.
BMC Systems Biology 4 (2010), 4.
}

\arguments{
  \item{links}{N x 2 matrix of link endpoints (N = number of links).}
  \item{samples}{Posterior samples, as given by either ICMg.combined.sampler or ICMg.links.sampler.}
}

\value{
  A matrix containing the component memberships for each data point (node).
}

\seealso{\code{\link{ICMg.combined.sampler}}, \code{\link{ICMg.links.sampler}}}
\keyword{ methods }

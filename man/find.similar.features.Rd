\name{find.similar.features}
\alias{find.similar.features}
\title{Find similar features with a given subnetwork.}
\description{Given subnetwork, orders the remaining features (genes)
in the input data based on similarity with the subnetwork. Allows the
identification of similar features that are not directly connected in
the input network.}
\usage{
find.similar.features(model, subnet.id, datamatrix = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{NetResponseModel object.}
  \item{subnet.id}{Investigated subnetwork.}
  \item{datamatrix}{Optional. Can be used to compare subnetwork
similarity with new data which was not used for learning the
subnetworks.}
  \item{verbose}{Logical indicating whether progress of the algorithm
should be indicated on the screen.}
}
\details{
The same similarity measure is used as when agglomerating the
subnetworks: the features are ordered by delta (change) in the cost
function, assuming that the feature would be merged in the subnetwork.
The smaller the change, the more similar the feature is (change would minimize the
new cost function value). Negative values of delta mean that the cost
function would be improved by merging the new feature in the
subnetwork, indicating features having coordinated response.
}
\value{A data frame with elements feature.names (e.g. gene IDs) and
delta, which indicates similarity level. See details for details. The
smaller, the more similar. The data frame is ordered such that the
features are listed by decreasing similarity.}
\references{ See citation("netresponse") for reference details.}
\author{ Leo Lahti \email{leo.lahti@iki.fi}}
%\note{}
%\seealso{}
\examples{
data(toydata)
model <- toydata$model
subnet.id <- "Subnet-1" 
g <- find.similar.features(model, subnet.id)
# List features that are similar to this subnetwork (delta < 0)
# (ordered by decreasing similarity)
subset(g, delta < 0)
}
\keyword{ utilities }


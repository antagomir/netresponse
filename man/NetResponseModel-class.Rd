\name{NetResponseModel-class}
%\Rdversion{1.1}
\docType{class}
\alias{NetResponseModel-class}

\alias{[[,NetResponseModel-method}
\alias{show,NetResponseModel-method}
%\alias{get.subnets}
%\alias{get.subnets,NetResponseModel-method} % NOTE: no space after comma!!

%\alias{get.model.parameters}
%\alias{get.model.parameters, NetResponseModel-method}

\title{Class "NetResponseModel"}
\description{A NetResponse model.}
\section{Objects from the Class}{
  Returned by \code{\link{detect.responses}} function.
}
\section{Slots}{
  \describe{
    \item{moves}{Subnetwork merging history.}
    %\item{costs}{Cost function values at the different steps.}
    \item{last.grouping}{Subnetworks in the last agglomeration level: feature indices}
    \item{subnets}{Subnetworks in the last agglomeration level: feature names}
    \item{params}{Input parameters.}
    \item{nodes}{Node ids from the original modelled data matrix.}
    \item{samples}{Sample ids from the original modelled data matrix. }
    \item{datamatrix}{Original input datamatrix that was used to learn the model.}
    \item{network}{Original input network that was used to learn the model.}
    \item{models}{Parameters for the learned subnetwork models.}
  }
}

%\section{Extends}{
%Class \code{"\linkS4class{list}"}, from data part.
%Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
%}

\section{Methods}{
  \describe{
    %\item{[}{\code{signature(x = "rpa")}: ... }
    \item{[[}{\code{signature(x = "NetResponseModel")}: ... }
    \item{show}{\code{signature(x = "NetResponseModel")}: ... }
%    \item{get.subnets}{\code{signature(model = "NetResponseModel",
%    get.names = "logical", stat = "data.frame",min.size =
%"numeric", max.size = "numeric", min.responses = "numeric")}: Returns
%a list of subnetworks detected by the model.}
}}


\author{Leo Lahti \email{leo.lahti@iki.fi}}

\examples{showClass("NetResponseModel")}

%\usage{\S4method{ get.subnets }{NetResponseModel}(model, level, get.names, stat, min.size, max.size, min.responses)}
%\seealso{\linkS4class{}}

%\examples{
%  #List the detected subnetworks (each is a list of nodes that belong to
%  #that subnetwork):
%  #res <- detect.responses(D, netw, verbose = FALSE)
%  #get.subnets(res)
%}
\keyword{classes}


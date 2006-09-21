\title{modifyWeights}
\name{modifyWeights}
\alias{modifyWeights}
\description{
Modify weights matrix for given gene status values.
}
\usage{
modifyWeights(weights=rep(1,length(status)), status, values, multipliers)
}
\arguments{
  \item{weights}{numeric matrix of relative weights, rows corresponding to genes and columns to arrays}
  \item{status}{character vector giving the control status of each spot on the array, of same length as the number of rows of \code{weights}}
  \item{values}{character vector giving subset of the unique values of \code{status}}
  \item{multipliers}{numeric vector of same length as \code{values} giving factor by which weights will be modified}
}

\details{
The function is usually used to temporarily modify the weights matrix during normalization of data.
The function can be used for example to give zero weight to spike-in ratio control spots during normalization.
}

\value{
Numeric matrix of same dimensions as \code{weights} with rows corresponding to \code{values} in \code{status} modified by the specified multipliers.
}
\author{Gordon Smyth}
\examples{
w <- matrix(runif(6*3),6,3)
status <- c("Gene","Gene","Ratio_Control","Ratio_Control","Gene","Gene")
modifyWeights(w,status,values="Ratio_Control",multipliers=0)
}
\seealso{
An overview of normalization functions available in LIMMA is given in \link{05.Normalization}.
}
\keyword{hplot}

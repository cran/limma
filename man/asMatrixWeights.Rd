\title{asMatrixWeights}
\name{asMatrixWeights}
\alias{asMatrixWeights}
\description{
Convert probe-weights or array-weights to a matrix of weights.
}
\usage{
asMatrixWeights(weights, dim)
}
\arguments{
  \item{weights}{numeric matrix of weights, rows corresponding to probes and columns to arrays.
  Or vector of probe weights. Or vector of array weights.}
  \item{dim}{numeric dimension vector of length 2, i.e., the number of probes and the number of arrays.}
}

\details{
This function converts a vector or probe-weights or a vector of array-weights to a matrix of the correct size.
Probe-weights are repeated across rows while array-weights are repeated down the columns.
If \code{weights} has length equal to the number of probes, it is assumed to contain probe-weights.
If it has length equal to the number of arrays, it is assumed to contain array-weights.
If the number of probes is equal to the number of arrays, then \code{weights} is assumed to contain array-weights if it is a row-vector of the correct size, i.e., if it is a matrix with one row.

This function is used internally by the linear model fitting functions in limma.
}

\value{
Numeric matrix of dimension \code{dim}.
}
\author{Gordon Smyth}
\examples{
asMatrixWeights(1:3,c(4,3))
asMatrixWeights(1:4,c(4,3))
}
\seealso{
\code{\link{modifyWeights}}.

An overview of functions in LIMMA used for fitting linear models is given in \link{06.LinearModels}.
}
\keyword{hplot}

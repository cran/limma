\name{exprs.MA}
\alias{exprs.MA}
\title{Extract Log-Expression Matrix from MAList}
\description{
Extract the matrix of log-expression values from an \code{MAList} object.
}
\usage{
exprs.MA(MA)
}
\arguments{
  \item{MA}{an \code{MAList} object.}
}
\details{
Converts M and A-values to log-expression values.
The output matrix will have two columns for each array, in the order green, red for each array.

This contrasts with \code{as.matrix.MAList} which extracts the M-values only, or \code{RG.MA} which converts to expression values in \code{RGList} form.
}
\value{
A numeric matrix with twice the columns of the input.
}
\author{Gordon Smyth}
\seealso{
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}
\keyword{array}

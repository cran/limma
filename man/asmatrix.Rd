\name{as.matrix}
\alias{as.matrix.RGList}
\alias{as.matrix.MAList}
\alias{as.matrix.MArrayLM}
\alias{as.matrix.marrayNorm}
\alias{as.matrix.exprSet}
\alias{as.matrix.PLMset}
\alias{as.matrix.ExpressionSet}
\alias{as.matrix.LumiBatch}
\alias{as.matrix.vsn}
\title{Turn a Microarray Data Object into a Matrix}
\description{
Turn a microarray data object into a numeric matrix by extracting the expression values.
}
\usage{
\method{as.matrix}{MAList}(x,\dots)
}
\arguments{
  \item{x}{an object of class \code{RGList}, \code{MAList}, \code{MArrayLM}, \code{marrayNorm}, \code{exprSet}, \code{PLMset}, \code{ExpressionSet}, \code{LumiBatch} or \code{vsn}.}
  \item{\dots}{additional arguments, not used for these methods.}
}
\details{
These methods extract the matrix of log-ratios, for \code{MAList} or \code{marrayNorm} objects, or the matrix of expression values for other expression objects such as \code{exprSet} or \code{ExressionSet}.
For \code{MArrayLM} objects, the matrix of fitted coefficients is extracted.

These methods involve loss of information, so the original data object is not recoverable.
}
\value{
A numeric matrix.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{as.matrix}} in the base package or \code{\link[Biobase]{exprs}} in the Biobase package.
  
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}
\keyword{array}

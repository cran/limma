\name{arrayWeights}
\alias{arrayWeights}
\title{Array Quality Weights}
\description{
Estimates relative quality weights for each array a multi-array experiment with replication.
}
\usage{
arrayWeights(y, fit)
}
\arguments{
  \item{y}{the data object used to estimate \code{fit}.
  Can be of any class which can be coerced to matrix, including \code{matrix}, \code{MAList}, \code{marrayNorm} or \code{exprSet}.}
  \item{fit}{\code{MArrayLM} fitted model object}
}
\details{
Estimates the relative reliability of each array by measuring how well the expression values for that array follow the linear model.
}
\value{
Numeric vector of weights, so same length as \code{ncol(fit)}.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{3.ReadingData}.
}
\examples{
\dontrun{
fit <- lmFit(y, design)
arrayWeights(y, fit)
}
}
\keyword{regression}

\name{ma3x3}
\alias{ma3x3.matrix}
\alias{ma3x3.spottedarray}
\title{Two dimensional Moving Averages with 3x3 Window}
\description{
Apply a specified function to each to each value of a matrix and its immediate neighbors.
}
\usage{
ma3x3.matrix(x,FUN=mean,na.rm=TRUE,...)
ma3x3.spottedarray(x,printer,FUN=mean,na.rm=TRUE,...)
}
\arguments{
  \item{x}{numeric matrix}
  \item{FUN}{function to apply to each window of values}
  \item{na.rm}{logical value, should missing values be removed when applying \code{FUN}}
  \item{...}{other arguments are passed to \code{FUN}}
  \item{printer}{list giving the printer layout, see \code{\link{PrintLayout-class}}}
}
\details{
For \code{ma3x3.matrix}, \code{x} is an arbitrary function.
for \code{ma3x3.spotted}, each column of \code{x} is assumed to contain the expression values of a spotted array in standard order.
The printer layout information is used to re-arrange the values of each column as a spatial matrix before applying \code{ma3x3.matrix}.
}
\value{
Numeric matrix of same dimension as \code{x} containing smoothed values
}
\seealso{
An overview of functions for background correction are given in \code{\link{04.Background}}.
}
\examples{
x <- matrix(c(2,5,3,1,6,3,10,12,4,6,4,8,2,1,9,0),4,4)
ma3x3.matrix(x,FUN="mean")
ma3x3.matrix(x,FUN="min")
}
\author{Gordon Smyth}
\keyword{smooth}

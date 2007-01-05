\name{avedups}
\alias{avedups}
\alias{avedups.default}
\alias{avedups.MAList}
\title{Average Over Duplicate Spots}
\description{
Condense a microarray data object so that values for within-array replicate spots are replaced with their average.
}
\usage{
\method{avedups}{default}(x, ndups=2, spacing=1, weights=NULL)
\method{avedups}{MAList}(x, ndups=x$printer$ndups, spacing=x$printer$spacing, weights=x$weights)
}
\arguments{
  \item{x}{a matrix-like object, usually a matrix or an \code{MAList} object.}
  \item{ndups}{number of within-array replicates for each probe.}
  \item{spacing}{number of spots to step from a probe to its duplicate.}
  \item{weights}{numeric matrix of spot weights.}
}
\details{
A new data object is computed in which each probe is represented by the (weighted) average of its duplicate spots.
For a \code{MAList} object, the components \code{M} and \code{A} are both averaged in this way.

If \code{x} is of mode \code{"character"}, then the duplicate values are assumed to be equal and the first is taken as the average.
}
\value{
A data object of the same class as \code{x} with \code{1/ndups} as many rows.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{dimnames}} in the base package.
  
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}
\keyword{array}

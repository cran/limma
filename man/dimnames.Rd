\name{dimnames}
\alias{dimnames.RGList}
\alias{dimnames.MAList}
\alias{dimnames.MArrayLM}
\title{Retrieve the Dimension Names of an RGList, MAList or MArrayLM Object}
\description{
Retrieve the dimension names of a microarray data object.
}
\usage{
\method{dimnames}{RGList}(x)
}
\arguments{
  \item{x}{an object of class \code{RGList}, \code{MAList} or \code{MArrayLM}}
}
\details{
The dimensionnames of an microarray object are the same as those of the most important matrix component of that object.

A consequece is that \code{rownames} and \code{colnames} will work as expected.
}
\value{
Either \code{NULL} or a list of length 2.
If a list, its components are either 'NULL' or a character vector the length of the appropriate dimension of \code{x}.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{dimnames}} in the base package.
  
  \link{02.Classes} gives an overview of data classes used in LIMMA.
}
\keyword{array}

\name{blockDiag}
\alias{blockDiag}

\title{Block Diagonal Matrix}
\description{Form a block diagonal matrix from the given blocks.}

\usage{
blockDiag(...)
}

\arguments{
\item{...}{numeric matrices}
}

\value{A block diagonal matrix with dimensions equal to the sum of the input dimensions}

\examples{
a <- matrix(1,3,2)
b <- matrix(2,2,2)
blockDiag(a,b)
}

\author{Gordon Smyth}

\seealso{
\link{10.Other}
}

\keyword{array}

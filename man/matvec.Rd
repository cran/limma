\name{matvec}
\alias{matvec}
\alias{vecmat}

\title{Multiply a Matrix by a Vector}
\description{Multiply the rows or columns of a matrix by the elements of a vector.}

\usage{
matvec(M, v)
vecmat(v, M)
}

\arguments{
\item{M}{numeric matrix, or object which can be coerced to a matrix.}
\item{v}{numeric vector, or object which can be coerced to a vector. Length should match the number of columns of \code{M} (for \code{matvec}) or the number of rows of \code{M} (for \code{vecmat})}
}

\value{A matrix of the same dimensions as \code{M}.}

\details{
\code{matvec(M,v)} is equivalent to \code{M \%*\% diag(v)} but is faster to execute.
Similarly \code{vecmat(v,M)} is equivalent to \code{diag(v) \%*\% M} but is faster to execute.
}

\examples{
A <- matrix(1:12,3,4)
A
matvec(A,c(1,2,3,4))
vecmat(c(1,2,3),A)
}

\author{Gordon Smyth}

\keyword{array}
\keyword{algebra}

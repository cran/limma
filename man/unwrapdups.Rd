\name{unwrapdups}
\alias{unwrapdups}
\title{Unwrap Duplicate Spot Values from Rows into Columns}
\description{Reshape a matrix so that a set of consecutive rows becomes a single row in the output.}
\usage{
unwrapdups(M,ndups=2,spacing=1)
}
\arguments{
  \item{M}{a matrix.}
  \item{ndups}{number of duplicate spots. The number of rows of M must be divisible by \code{ndups}.}
  \item{spacing}{the spacing between the rows of \code{M} corresponding to duplicate spots, \code{spacing=1} for consecutive spots}
}
\value{A matrix containing the same values as \code{M} but with fewer rows and more columns by a factor of \code{ndups}.
Each set of \code{ndups} rows in \code{M} is strung out to a single row so that duplicate values originally in consecutive rows in the same column are in consecutive columns in the output.
}
\details{
This function is used on matrices corresponding to a series of microarray experiments.
Rows corresponding to duplicate spots are re-arranged to that all values corresponding to a single gene are on the same row.
This facilitates fitting models or computing statistics for each gene.
}
\author{Gordon Smyth}
\examples{
M <- matrix(1:12,6,2)
unwrapdups(M,ndups=2)
unwrapdups(M,ndups=3)
unwrapdups(M,ndups=2,spacing=3)
}
\keyword{array}

\name{is.fullrank}
\alias{is.fullrank}
\alias{nonEstimable}

\title{Check for Full Column Rank}

\description{
Test whether a numeric matrix has full column rank.
}

\usage{
is.fullrank(x)
nonEstimable(x)
}

\arguments{
\item{x}{a numeric matrix or vector}
}

\value{
\code{is.fullrank} returns \code{TRUE} or \code{FALSE}.

\code{nonEstimable} returns a character vector of names for the columns of \code{x} which are linearly dependent on previous columns.
If \code{x} has full column rank, then the value is \code{NULL}.
}

\details{
\code{is.fullrank} is used to check the integrity of design matrices in limma, for example after \link[limma:subsetting]{subsetting} operations.

\code{nonEstimable} is used by \code{\link{lmFit}} to report which coefficients in a linear model cannot be estimated.
}

\author{Gordon Smyth}

\examples{
# TRUE
is.fullrank(1)
is.fullrank(cbind(1,0:1))

# FALSE
is.fullrank(0)
is.fullrank(matrix(1,2,2))
nonEstimable(matrix(1,2,2))
}
\keyword{algebra}

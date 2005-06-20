\name{printHead}
\alias{printHead}
\title{Print Leading Rows of Large Objects}

\description{
Print the leading rows of a large vector, matrix or data.frame.
This function is used by \code{show} methods for data classes defined in LIMMA.
}

\usage{
printHead(x)
}

\arguments{
  \item{x}{any object}
}

\details{
If \code{x} is a vector with more than 20 elements, then \code{printHead(x)} prints only the first 5 elements.
If \code{x} is a matrix or data.frame with more than 10 rows, then \code{printHead(x)} prints only the first 5 rows.
Any other type of object is printed normally.
}

\author{Gordon Smyth}

\seealso{
An overview of classes defined in LIMMA is given in \link{02.Classes}
}

\keyword{hplot}

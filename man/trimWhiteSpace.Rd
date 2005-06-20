\name{trimWhiteSpace}
\alias{trimWhiteSpace}

\title{Trim Leading and Trailing White Space}
\description{Trims leading and trailing white space from character strings.}

\usage{
trimWhiteSpace(x)
}

\arguments{
\item{x}{character vector}
}

\value{
A character vector of the same length as \code{x} in which leading and trailing white space has been stripped off each value.
}

\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\examples{
x <- c("a "," b ")
trimWhiteSpace(x)
}

\author{Tim Beissbarth and Gordon Smyth}

\keyword{character}

\name{isNumeric}
\alias{isNumeric}

\title{Test for Numeric Argument}
\description{
Test whether argument is numeric or a data.frame with numeric columns.
}

\usage{
isNumeric(x)
}

\arguments{
\item{x}{any object}
}

\value{\code{TRUE} or \code{FALSE}}

\details{
This function is used to check the validity of arguments for numeric functions.
It is an attempt to emulate the behavior of internal generic math functions.

\code{isNumeric} differs from \code{is.numeric} in that data.frames with all columns numeric are accepted as numeric.
}

\author{Gordon Smyth}

\examples{
isNumeric(3)
isNumeric("a")
x <- data.frame(a=c(1,1),b=c(0,1))
isNumeric(x)   # TRUE
is.numeric(x)  # FALSE
}
\seealso{
   \code{\link[base]{is.numeric}}, \code{\link[base]{Math}}
}
\keyword{programming}

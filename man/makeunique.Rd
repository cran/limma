\name{makeUnique}
\alias{makeUnique}
\title{Make Values of Character Vector Unique}
\description{
Paste characters on to values of a character vector to make them unique.
}
\usage{
makeUnique(x)
}
\arguments{
  \item{x}{object to be coerced to a character vector}
}
\details{
Repeat values of \code{x} are labelled with suffixes "1", "2" etc.
}
\value{
A character vector of the same length as \code{x}
}
\author{Gordon Smyth}
\seealso{
\code{makeUnique} is called by \code{\link{merge.RGList}}.
Compare with \code{\link[base]{make.unique}} in the base package.
}
\examples{
x <- c("a","a","b")
makeUnique(x)
}
\keyword{character}

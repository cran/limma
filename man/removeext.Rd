\name{removeExt}
\alias{removeExt}

\title{Remove Common Extension from File Names}
\description{Finds and removes any common extension from a vector of file names.}

\usage{
removeExt(x)
}

\arguments{
\item{x}{character vector}
}

\value{
A character vector of the same length as \code{x} in which any common extension has been stripped off.
}

\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\examples{
x <- c("slide1.spot","slide2.spot","slide3.spot")
removeExt(x)
}

\author{Gordon Smyth}

\keyword{character}
\keyword{file}

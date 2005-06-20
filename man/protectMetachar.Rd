\name{protectMetachar}
\alias{protectMetachar}

\title{Protect Metacharacters}
\description{Add backslashes before any metacharacters found in a string.}

\usage{
protectMetachar(x)
}

\arguments{
\item{x}{character vector}
}

\details{
This function is used to protect strings containing metacharacters so that the metacharacters can be treated as ordinary characters in string matching functions operations.
}

\value{
A character vector of the same length as \code{x} in which two backslashes have been inserted before any metacharacter.
}

\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\examples{
# without protectMetachar, this would be no match
grep(protectMetachar("Ch1 (mean)"),"Ch1 (mean)")
}

\author{Gordon Smyth}

\keyword{character}

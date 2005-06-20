\name{LargeDataObject-class}
\docType{class}
\alias{LargeDataObject-class}
\alias{show,LargeDataObject-method}
\title{Large Data Object - class}

\description{
A virtual class including the data classes \code{RGList}, \code{MAList} and \code{MArrayLM}, all of which typically contain large quantities of numerical data in vector, matrices and data.frames.
}

\section{Methods}{
A \code{show} method is defined for objects of class \code{LargeDataObject} which uses \code{printHead} to print only the leading elements or rows of components or slots which contain large quantities of data.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
}

\examples{
#  see normalizeBetweenArrays
}

\keyword{classes}
\keyword{data}

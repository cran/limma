\name{getSpacing}
\alias{getSpacing}
\title{Get Numerical Spacing}
\description{
Convert character to numerical spacing measure for within-array replicate spots.
}
\usage{
getSpacing(spacing, layout)
}
\arguments{
  \item{spacing}{character string or integer.
  Acceptable character strings are \code{columns}, \code{rows} or \code{topbottom}.
  Integer values are simply passed through.}
  \item{layout}{list containing printer layout information}
}
\value{
Integer giving spacing between replicate spots in the gene list.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\examples{
getSpacing("columns",list(ngrid.r=2,ngrid.c=2,nspot.r=20,nspot.c=19))
getSpacing("rows",list(ngrid.r=2,ngrid.c=2,nspot.r=20,nspot.c=19))
getSpacing("topbottom",list(ngrid.r=2,ngrid.c=2,nspot.r=20,nspot.c=19))
}
\keyword{IO}

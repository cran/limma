\name{PrintLayout}
\docType{class}
\alias{PrintLayout-class}
\title{Print Layout - class}

\description{
A list-based class for storing information about the process used to print spots on a microarray.

\code{PrintLayout} objects can be created using \code{\link{getLayout}}.
The \code{printer} component of an \code{RGList} or \code{MAList} object is of this class.
}

\section{Slots/List Components}{
Objects of this class contains no slots but should contain the following list components:
\tabular{ll}{
  \code{ngrid.r}:\tab number of grid rows on the arrays\cr
  \code{ngrid.c}:\tab number of grid columns on the arrays\cr
  \code{nspot.r}:\tab number of rows of spots in each grid\cr
  \code{nspot.c}:\tab number of columns of spots in each grid\cr
  \code{ndups}:\tab number of duplicates of each DNA clone, i.e., number of times print-head dips into each well of DNA\cr
  \code{spacing}:\tab number of spots between duplicate spots.  Only applicable if \code{ndups>1}.
  \code{spacing=1} for side-by-side spots by rows, \code{spacing=nspot.c} for side-by-side spots by columns, \code{spacing=ngrid.r*ngrid.c*nspot.r*nspot.c/2} for duplicate spots in top and bottom halves of each array.\cr
  \code{npins}:\tab actual number of pins or tips on the print-head\cr
  \code{start}:\tab character string giving position of the spot printed first in each grid.
  Choices are \code{"topleft"} or \code{"topright"} and partial matches are accepted.
}
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
}

\examples{
#  Settings for Swirl and ApoAI example data sets in User's Guide

printer <- list(ngrid.r=4, ngrid.c=4, nspot.r=22, nspot.c=24, ndups=1, spacing=1, npins=16, start="topleft")

#  Typical settings at the Australian Genome Research Facility

#  Full pin set, duplicates side-by-side on same row
printer <- list(ngrid.r=12, ngrid.c=4, nspot.r=20, nspot.c=20, ndups=2, spacing=1, npins=48, start="topright")

#  Half pin set, duplicates in top and lower half of slide
printer <- list(ngrid.r=12, ngrid.c=4, nspot.r=20, nspot.c=20, ndups=2, spacing=9600, npins=24, start="topright")
}

\keyword{classes}
\keyword{data}

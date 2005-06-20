\name{printorder}
\alias{printorder}
\title{Identify Order in which Spots were Printed}
\description{
Identify order in which spots were printed and the 384-well plate from which they were printed.
}
\usage{
printorder(layout, ndups=1, spacing="columns", npins, start="topleft")
}
\arguments{
  \item{layout}{list with the components \code{ngrid.r}, \code{ngrid.c}, \code{nspot.r} and \code{nspot.c}, or an \code{RGList} or \code{MAList} object from which the printer layout may be extracted.}
  \item{ndups}{number of duplicate spots, i.e., number of times print-head dips into each well}
  \item{spacing}{character string indicating layout of duplicate spots.
  Choices are \code{"columns"}, \code{"rows"} or \code{"topbottom"}.}
  \item{npins}{actual number of pins or tips on the print-head}
  \item{start}{character string giving position of the spot printed first in each grid.
  Choices are \code{"topleft"} or \code{"topright"} and partial matches are accepted.}
}
\details{
In most cases the printer-head contains the \code{layout$ngrid.r} times \code{layout$ngrid.c} pins or tips and the array is printed using \code{layout$nspot.r} times \code{layout$npot.c} dips of the head.
The plate holding the DNA to be printed is assumed to have 384 wells in 16 rows and 24 columns.

\code{ndups} indicates the number of spots printed from each well.
The replicate spots from multiple dips into the same wells are assumed to be side-by-side by columns (\code{spacing="columns"}), by rows (\code{spacing="rows"}) or in the top and bottom halves of the array (\code{spacing="topbottom"}).

In some cases a smaller number of physical pins is used and the total number of grids is built up by effectively printing two or more sub-arrays on the same slide.
In this case the number of grids should be a multiple of the number of pins.

Printing is assumed to proceed by rows within in each grid starting either from the top-left or the top-right.
}
\value{
List with components
\item{printorder}{numeric vector giving printorder of each spot, i.e., which dip of the print-head was used to print it}
\item{plate}{numeric vector giving plate number from which each spot was printed}
\item{plate.r}{numeric vector giving plate-row number of the well from which each spot was printed}
\item{plate.c}{numeric vector giving plate-column number of the well from which each spot was printed}
\item{plateposition}{character vector summarizing plate number and plate position of the well from which each spot was printed with letters for plate rows and number for columns.
For example \code{02B13} is second row, 13th column, of the second plate.}
}
\seealso{
\code{\link{normalizeForPrintorder}}.

An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\author{Gordon Smyth}
\examples{
printorder(list(ngrid.r=2,ngrid.c=2,nspot.r=12,nspot.c=8))
}
\keyword{IO}

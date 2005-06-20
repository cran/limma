\name{getLayout}
\alias{getLayout}
\alias{getLayout2}
\title{Extract the Print Layout of an Array from the GAL File}
\description{
From the Block, Row and Column information in a genelist, determine the number of grid rows and columns on the array and the number of spot rows and columns within each grid.
}
\usage{
getLayout(gal, guessdups=FALSE)
getLayout2(galfile)
}
\arguments{
  \item{gal}{data.frame containing the GAL, i.e., giving the position and gene identifier of each spot}
  \item{galfile}{name or path of GAL file}
  \item{guessdups}{logical, if \code{TRUE} then try to determine number and spacing of duplicate spots, i.e., within-array replicates}
}
\details{
A GenePix Array List (GAL) file is a list of genes and associated information produced by an Axon microarray scanner.
The function \code{getLayout} determines the print layout from a data frame created from a GAL file or gene list.
The data.frame must contain columns \code{Block}, \code{Column} and \code{Row}.
(The number of tip columns is not assumed to be either one or four.)
The function \code{getLayout2} determines the print layout from the header information of an actual GAL file.
}
\value{
A \code{printlayout} object, which is a list with the following components.
The last two components are present only if \code{guessdups=TRUE}.
  \item{ngrid.r}{integer, number of grid rows on the arrays}
  \item{ngrid.c}{integer, number of grid columns on the arrays}
  \item{nspot.r}{integer, number of rows of spots in each grid}
  \item{nspot.c}{integer, number of columns of spots in each grid}
  \item{ndups}{integer, number of times each probe is printed on the array}
  \item{spacing}{integer, spacing between multiple printings of each probe}
}
\author{Gordon Smyth and James Wettenhall}
\seealso{
\code{\link[marray:marrayTools]{gpTools}}.

An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\examples{
# gal <- readGAL()
# layout <- getLayout(gal)
}
\keyword{IO}

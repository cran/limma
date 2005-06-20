\name{readGAL}
\alias{readGAL}
\title{Read a GAL file}
\description{
Read a GenePix Array List (GAL) file into a dataframe.
}
\usage{
readGAL(galfile=NULL,path=NULL,header=TRUE,sep="\t",quote="\"",skip=NULL,as.is=TRUE,...)
}
\arguments{
  \item{galfile}{character string giving the name of the GAL file.  If \code{NULL} then a file with extension \code{.gal} is found in the directory specified by \code{path}.}
  \item{path}{character string giving the directory containing the files.  If \code{NULL} then assumed to be the current working directory.}
  \item{header}{logical variable, if \code{TRUE} then the first line after \code{skip} is assumed to contain column headings.  If \code{FALSE} then a value should specified for \code{skip}.}
  \item{sep}{the field separator character}
  \item{quote}{the set of quoting characters}
  \item{skip}{number of lines of the GAL file to skip before reading data.  If \code{NULL} then this number is determined by searching the file for column headings.}
  \item{as.is}{logical variable, if \code{TRUE} then read in character columns as vectors rather than factors.}
  \item{...}{any other arguments are passed to \code{read.table}}
}
\details{
A GAL file is a list of genes IDs and associated information produced by an Axon microarray scanner.
Apart from header information, the file must contain data columns labeled \code{Block}, \code{Column}, \code{Row} and \code{ID}.
A \code{Name} column is usually included as well.
Other columns are optional.
See the Axon URL below for a detaile description of the GAL file format.

This function reads in the data columns with a minimum of user information.
In most cases the function can be used without specifying any of the arguments.
}
\value{
A data frame with columns
  \item{Block}{numeric vector containing the print tip indices}
  \item{Column}{numeric vector containing the spot columns}
  \item{Row}{numeric vector containing the spot rows}
  \item{ID}{character vector, for factor if \code{as.is=FALSE}, containing gene library identifiers}
  \item{Name}{character vector, for factor if \code{as.is=FALSE}, containing gene names}
The data frame will be sorted so that \code{Column} is the fastest moving index, then \code{Row}, then \code{Block}.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\references{
\url{http://www.axon.com/gn_GenePix_File_Formats.html}
}
\examples{
# readGAL()
# will read in the first GAL file (with suffix ".gal")
# found in the current working directory
}
\keyword{IO}

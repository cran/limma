\name{readImaGeneHeader}
\alias{readImaGeneHeader}
\title{Read ImaGene Header Information}
\description{
Read the header information from an ImaGene image analysis output file.
This function is used internally by \code{read.maimages} and is not usually called directly by users.
}
\usage{
readImaGeneHeader(file)
}
\arguments{
  \item{file}{character string giving file name or path}
}
\details{
The raw data files exported by the image analysis software ImaGene include a number of header lines which contain information about the printing process.
This function extracts that information and locates the line where the intensity data begins.
}
\value{
A list with components
  \item{NHeaderRecords}{number of lines in the file before the intensity data begins}
  \item{Version}{version number of ImaGene software}
  \item{Date}{character string giving time and data that array was scanned}
  \item{ImageFile}{character string giving original file name to which data was written by ImaGene}
  \item{Inverted}{logical}
  \item{FieldDimension}{list with components
    \code{Field} containing a character string,
    \code{Metarows} containing number of grid rows,
    \code{Metacols} containing number of grid columns,
    \code{Rows} containing number of spot rows in each grid,
    \code{Cols} containing number of spot columns in each grid}
  \item{MeasurementParameters}{list with numerical components \code{SignalLow}, \code{SignalHigh}, \code{BackgroundLow}, \code{BackgroundHigh}, \code{BackgroundBuffer} and \code{BackgroundWidth}}
}
\references{
\url{http://www.biodiscovery.com/imagene.asp}
}
\author{Gordon Smyth}
\seealso{\code{\link{read.imagene}}

An overview of LIMMA functions to read data is given in \link{03.ReadingData}.
}
\keyword{file}

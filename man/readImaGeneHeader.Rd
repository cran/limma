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
The raw data files exported by the image analysis software ImaGene include a number of header lines which contain information about the printing and scanning processes.
This function extracts that information and locates the line where the intensity data begins.
}
\value{
A list containing information read from the header of the ImaGene file.
Each Begin-End environment found in the file header will become a recursive list in the output object, with components corresponding to fields in the file.
See the ImaGene documentation for further information.
The output object will also contain a component \code{NHeaderRecords}giving the number of lines in the file before the intensity data begins.
}
\references{
\url{http://www.biodiscovery.com/imagene.asp}
}
\author{Gordon Smyth}
\seealso{\code{\link{read.imagene}}

An overview of LIMMA functions to read data is given in \link{03.ReadingData}.
}
\examples{\dontrun{
h <- readImaGeneHeader("myImaGenefile.txt")
names(h)
h$NHeaderRecords
h[["Field Dimensions"]]
}
}
\keyword{file}

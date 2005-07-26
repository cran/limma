\name{readHeader}
\alias{readBlueFuseHeader}
\alias{readGPRHeader}
\alias{readSMDHeader}
\title{Read Header Information from Image Analysis Raw Data File}
\description{
Read the header information from a GenePix Results (GPR) file or from an SMD raw data file.
These functions are used internally by \code{read.maimages} and are not usually called directly by users.
}
\usage{
readBlueFuseHeader(file)
readGPRHeader(file)
readSMDHeader(file)
}
\arguments{
  \item{file}{character string giving file name. If it does not contain an absolute path, the file name is relative to the current working directory.}
}
\details{
Raw data files exported by image analysis programs include a number of header lines which contain information about the scanning process.
This function extracts that information and locates the line where the intensity data begins.
\code{readBlueFuseHeader} is for BlueFuse output, \code{readGPRHeader} is for GenePix output, \code{readSMDHeader} is for files from the Stanford Microarray Database (SMD).
}
\value{
A list with components corresponds to lines of header information.
A key component is \code{NHeaderRecords} which gives the number of lines in the file before the intensity data begins.
All other components are character vectors.
}
\references{
See \url{http://www.axon.com/gn_GenePix_File_Formats.html} for GenePix formats.

See \url{http://www.bluegnome.co.uk} for information on BlueFuse.

See \url{http://genome-www.stanford.edu/Microarray} for the SMD.
}
\author{Gordon Smyth. Gregory Lefebvre wrote the original version of readBlueFuseHeader.}
\seealso{\code{\link{read.maimages}}

An overview of LIMMA functions to read data is given in \link{03.ReadingData}.
}
\keyword{file}

\name{readHeader}
\alias{readGPRHeader}
\alias{readSMDHeader}
\title{Read Header Information from Image Analysis Raw Data File}
\description{
Read the header information from a GenePix Results (GPR) file or from an SMD raw data file.
These functions are used internally by \code{read.maimages} and is not usually called directly by users.
}
\usage{
readGPRHeader(file)
readSMDHeader(file)
}
\arguments{
  \item{file}{character string giving file name or including path}
}
\details{
Raw data files exported by image analysis programs include a number of header lines which contain information about the scanning process.
This function extracts that information and locates the line where the intensity data begins.
\code{readGPRHeader} is for GenePix output, \code{readSMDHeader} is for 
}
\value{
A list with components corresponds to lines of header information.
All components are character vectors.
A key component is \code{NHeaderRecords} which gives the number of lines in the file before the intensity data begins.
}
\references{
\url{http://www.axon.com/gn_GenePix_File_Formats.html}

\url{http://genome-www.stanford.edu/Microarray}
}
\author{Gordon Smyth}
\seealso{\code{\link{read.maimages}}

An overview of LIMMA functions to read data is given in \link{03.ReadingData}.
}
\keyword{file}

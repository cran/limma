\name{read.matrix}
\alias{read.matrix}
\title{Read Matrix with Headers from File}
\description{
Read a numeric matrix from a file assuming column headings on the first line.
Not normally used directly by users.

This function is deprecated in limma 2.9.7 (1 January 2007).
}
\usage{
read.matrix(file,nrows=0,skip=0,...)
}
\arguments{
  \item{file}{character string giving the file name}
  \item{nrows}{maximum number of rows of data to read, if greater than zero}
  \item{skip}{number of lines of the data file to skip before reading data}
  \item{...}{any other arguments to be passed to \code{scan}}
}
\details{
This function is similar to but faster than \code{read.table(file,header=TRUE)} when all the columns are numeric.
}
\value{
A numeric matrix with column names.
}
\author{Gordon Smyth}
\seealso{
\code{\link[base]{read.table}}, \code{\link[base]{scan}}.

An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\keyword{file}

\name{readTargets}
\alias{readTargets}
\title{Read RNA Targets File}
\description{
Read targets file for a two-color microarray experiment into a dataframe.
}
\usage{
readTargets(file="Targets.txt", path=NULL, sep="\t", row.names="FileName")
}
\arguments{
  \item{file}{character string giving the name of the targets file.}
  \item{path}{character string giving the directory containing the file.
       Can be omitted if the file is in the current working irectory.} 
  \item{sep}{field separator character}
  \item{row.names}{character string giving name of column from which to obtain row names}
}
\details{
The targets file is a text file containing information about the RNA samples used as targets in the microarray experiment.
The files should have a row of column headers and shoule have one row for each array in the experiment.
Normally there are columns labelled \code{Cy3} and \code{Cy5} or similar specifying which RNA samples are hybridized to each channel of each array.
Other columns may contain any other information associated with the arrays or targets used in the experiment.
If there is a column named \code{FileName} then that will be used to construct row names for the dataframe.

See the Limma User's Guide for examples of this function.
}
\value{
Data frame usually including columns
  \item{Cy3}{character vector giving names of RNA samples}
  \item{Cy5}{character vector giving names of RNA samples}
Character columns are not converted into factors.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{3.ReadingData}.
}
\keyword{IO}

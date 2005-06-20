\name{readTargets}
\alias{readTargets}
\title{Read Targets File}
\description{
Read targets file for a microarray experiment into a dataframe.
}
\usage{
readTargets(file="Targets.txt", path=NULL, sep="\t", row.names=NULL, quote="\"",...)
}
\arguments{
  \item{file}{character string giving the name of the targets file.}
  \item{path}{character string giving the directory containing the file.
       Can be omitted if the file is in the current working irectory.} 
  \item{sep}{field separator character}
  \item{row.names}{character string giving the name of a column from which to obtain row names}
  \item{quote}{the set of quoting characters}
  \item{...}{other arguments are passed to \code{\link[base]{read.table}}}
}
\details{
The targets file is a text file containing information about the RNA samples used as targets in the microarray experiment.
Rows correspond to arrays and columns to covariates associated with the targets.
For a two-color experiment, the targets file will normally include columns labelled \code{Cy3} and \code{Cy5} or similar specifying which RNA samples are hybridized to each channel of each array.
Other columns may contain any other covariate information associated with the arrays or targets used in the experiment.

If \code{row.names} is non-null and there is a column by that name with unique values, then those values will be used as row names for the dataframe.
If \code{row.names} is null, then the column \code{Labels} will be used if such exists or, failing that, the column \code{FileName}.

See the Limma User's Guide for examples of this function.
}
\value{
A dataframe. Character columns are not converted into factors.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\keyword{IO}

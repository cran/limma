\name{readSpotTypes}
\alias{readSpotTypes}
\title{Read Spot Types File}
\description{
Read a table giving regular expressions to identify different types of spots in the gene-dataframe.
}
\usage{
readSpotTypes(file="SpotTypes.txt",path=NULL,sep="\t",check.names=FALSE,...)
}
\arguments{
  \item{file}{character string giving the name of the file specifying the spot types.}
  \item{path}{character string giving the directory containing the file.
       Can be omitted if the file is in the current working irectory.} 
  \item{sep}{the field separator character}
  \item{check.names}{logical, if \code{FALSE} column names will not be converted to valid variable names, for example spaces in column names will not be left as is}
  \item{\dots}{any other arguments are passed to \code{read.table}}
}
\details{
The file is a text file with rows corresponding to types of spots and the following columns: \code{SpotType} gives the name for the spot type, \code{ID} is a regular expression matching the ID column, \code{Name} is a regular expression matching the Name column, and \code{Color} is the R name for the color to be associated with this type.
}
\value{
A data frame with columns
  \item{SpotType}{character vector giving names of the spot types}
  \item{ID}{character vector giving regular expressions}
  \item{Name}{character vector giving regular expressions}
  \item{Color}{character vector giving names of colors}
}
\author{Gordon Smyth following idea of James Wettenhall}
\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\keyword{IO}

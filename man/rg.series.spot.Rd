\name{rg.series.spot}
\alias{rg.series.spot}
\title{Read RGList from SPOT Image Analysis Output Files}
\description{
Extracts an RGList from a series of Spot image analysis files.

This is a depreciated function.
Use \code{\link{read.maimages}} instead.
}
\usage{
rg.series.spot(slides,path=NULL,names.slides=names(slides),suffix="spot",wt.fun=NULL,verbose=TRUE,...)
}
\arguments{
  \item{slides}{character vector giving the names of the files containing the Spot output}
  \item{path}{character string giving the directory containing the files. Can be omitted if the files are in the current working directory.}
  \item{names.slides}{names to be associated with each slide as column name}
  \item{suffix}{the file names are assumed to have this suffix added to names in \code{slides}}
  \item{wt.fun}{function to calculate quality weights}
  \item{verbose}{\code{TRUE} to report each time a file in read in}
  \item{\dots}{any other arguments to be passed to \code{scan}}
}
\details{
This function extracts the foreground and background intensities from a series of files, produced by the image analysis program SPOT, and assembles them in the components of one list.

Spot quality weights may also be extracted using an optional weight function.
}
\value{
A list containing the components
  \item{R}{matrix containing the red channel foreground intensities for each spot for each array.}
  \item{Rb}{matrix containing the red channel background intensities for each spot for each array.}
  \item{G}{matrix containing the green channel foreground intensities for each spot for each array.}
  \item{Gb}{matrix containing the green channel background intensities for each spot for each array.}
  \item{weights}{spot quality weights, if \code{wt.fun} is not \code{NULL}}
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\keyword{IO}

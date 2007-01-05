\name{rg.spot}
\alias{rg.spot}
\title{Extract RGList from data.frames Containing SPOT Data}
\description{
Extracts an RGList from Spot image analysis output when the data has already been read from files into data.frames objects.

This function is deprecated in limma 2.9.7 (1 January 2007).
}
\usage{
rg.spot(slides,names.slides=names(slides),suffix="spot",area=FALSE)
}
\arguments{
  \item{slides}{Character vector giving the names of the data frames containing the Spot output.}
  \item{names.slides}{Names to be associated with each slide as column name.}
  \item{suffix}{The dataframe names are assumed to have this suffix added to names in \code{slides}.}
  \item{area}{If \code{TRUE} then the output list includes a component containing the spot areas.}
}
\details{
This function extracts the foreground and background intensities from a series of data frames and assembles them in the components of one list.

Spot areas may also be extracted, which is useful for downweighting unusually small or large spots in subsequent analyses.
}
\value{
A list containing the components
  \item{R}{A matrix containing the red channel foreground intensities for each spot for each array.}
  \item{Rb}{A matrix containing the red channel background intensities for each spot for each array.}
  \item{G}{A matrix containing the green channel foreground intensities for each spot for each array.}
  \item{Gb}{A matrix containing the green channel background intensities for each spot for each array.}
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\keyword{IO}

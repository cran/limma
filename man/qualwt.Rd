\name{QualityWeights}
\alias{QualityWeights}
\alias{wtarea}
\alias{wtflags}
\alias{wtIgnore.Filter}
\title{Spot Quality Weights}
\description{
Functions to calculate quality weights for individual spots based on image analyis output file.
}
\usage{
wtarea(ideal=c(160,170))
wtflags(weight=0,cutoff=0)
wtIgnore.Filter
}
\arguments{
  \item{ideal}{numeric vector giving the ideal area or range of areas for a spot in pixels}
  \item{weight}{weight to be given to flagged spots}
  \item{cutoff}{cutoff value for \code{Flags} below which spots will be downweighted}
}
\details{
These functions can be passed as an argument to \code{read.maimages} to construct quality weights as the microarray data is read in.

\code{wtarea} downweights unusually small or large spots and is designed for SPOT output.
It gives weight 1 to spots which have areas in the ideal range, given in pixels, and linearly downweights spots which are smaller or larger than this range.

\code{wtflags} is designed for GenePix output and gives the specified weight to spots with \code{Flags} value less than the \code{cutoff} value.
Choose \code{cutoff=0} to downweight all flagged spots.
Choose \code{cutoff=-50} to downweight bad or absent spots or \code{cutoff=-75} to downweight only spots which have been manually flagged as bad.

\code{wtIgnore.Filter} is designed for QuantArray output and sets the weights equal to the column \code{Ignore Filter} produced by QuantArray.
These weights are 0 for spots to be ignored and 1 otherwise.
}
\value{
A function which takes a dataframe or matrix as argument and produces a numeric vector of weights between 0 and 1
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\examples{
#  Read in spot output files from current directory and give full weight to 165
#  pixel spots.  Note: for this example to run you must set fnames to the names
#  of actual spot output files (data not provided).
\dontrun{
RG <- read.maimages(fnames,source="spot",wt.fun=wtarea(165))
#  Spot will be downweighted according to weights found in RG
MA <- normalizeWithinArrays(RG,layout)
}
}
\keyword{regression}

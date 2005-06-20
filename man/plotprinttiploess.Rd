\title{MA Plots by Print-Tip Group}
\name{plotPrintTipLoess}
\alias{plotPrintTipLoess}
\description{
Creates a coplot giving MA-plots with loess curves by print-tip groups.
}
\usage{
plotPrintTipLoess(object,layout,array=1,span=0.4,...)
}
\arguments{
  \item{object}{\code{MAList} or \code{RGList} object or list with components \code{M} containing log-ratios and \code{A} containing average intensities}
  \item{layout}{a list specifying the number of tip rows and columns and the number of spot rows and columns printed by each tip.
  Defaults to \code{MA$printer} if that is non-null.}
  \item{array}{integer giving the array to be plotted. Corresponds to columns of \code{M} and \code{A}.}
  \item{span}{span of window for \code{lowess} curve}
  \item{...}{other arguments passed to \code{panel.smooth}}
}
\details{
Note that spot quality weights in \code{object} are not used for computing the loess curves for this plot even though such weights would be used for loess normalization using \code{normalizeWithinArrays}.
}
\value{
A plot is created on the current graphics device.
If there are missing values in the data, then the vector of row numbers for spots with missing values is invisibly returned, as for \code{coplot}.
}
\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\author{Gordon Smyth}
\keyword{hplot}

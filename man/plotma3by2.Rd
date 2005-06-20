\title{Write MA-Plots to Files}
\name{plotMA3by2}
\alias{plotMA3by2}
\description{
Write MA-plots to files in PNG format, six plots to a file in a 3 by 2 grid arrangement.
}
\usage{
plotMA3by2(MA, prefix="MA", main=colnames(MA), zero.weights=FALSE, common.lim=TRUE, ...)
}
\arguments{
  \item{MA}{an \code{MAList} or \code{RGList} object, or any list with components \code{M} containing log-ratios and \code{A} containing average intensities}
  \item{prefix}{character string giving prefix to attach to file names}
  \item{main}{character vector giving titles for plots}
  \item{zero.weights}{logical, should points with non-positive weights be plotted}
  \item{common.lim}{logical, should all plots on a page use the same axis limits}
  \item{...}{any other arguments are passed to \code{plotMA}}
}

\details{
This function writes plots in PNG format in an arrangement optimized for A4-sized paper.
}

\value{
No value is returned, but one or more files are written to the working directory.
The number of files is determined by the number of columns of \code{MA}.
}

\author{Gordon Smyth}

\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\keyword{hplot}

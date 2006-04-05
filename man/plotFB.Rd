\title{FB-Plot}
\name{plotFB}
\alias{plotFB}
\description{
Creates foreground-background plots.
}
\usage{
plotFB(RG, array=1, lim="separate", pch=16, cex=0.2, ...)
}
\arguments{
  \item{RG}{an \code{RGList} object.}
  \item{array}{integer giving the array to be plotted. Corresponds to columns of \code{R}, \code{G}, \code{Rb} and \code{Gb}.}
  \item{lim}{character string indicating whether the red and green plots should have \code{"separate"} or \code{"common"} x- and y- co-ordinate limits.}
  \item{pch}{vector or list of plotting characters. Defaults to integer code 16.}
  \item{cex}{numeric vector of plot symbol expansions.} 
  \item{...}{any other arguments are passed to \code{plot}}
}

\details{
A foreground-background plot is a plot of log2-foreground vs log2-background for a particular channel on a particular two-color array.
This function produces a pair of plots, one for green and one for red, for a specified array.

See \code{\link[graphics]{points}} for possible values for \code{pch}, \code{col} and \code{cex}.
}

\value{A plot is created on the current graphics device.}
\author{Gordon Smyth}
\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\keyword{hplot}

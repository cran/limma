\title{Write Imageplots to Files}
\name{imageplot3by2}
\alias{imageplot3by2}
\description{
Write imageplots to files in PNG format, six plots to a file in a 3 by 2 grid arrangement.
}
\usage{
imageplot3by2(RG, z="Gb", prefix=paste("image",z,sep="-"), path=NULL, zlim=NULL, common.lim=TRUE, ...)
}
\arguments{
  \item{RG}{an \code{RGList} or \code{MAList} object, or any list with component named by \code{z}}
  \item{z}{character string giving name of component of \code{RG} to plot}
  \item{prefix}{character string giving prefix to attach to file names}
  \item{path}{character string specifying directory for output files}
  \item{zlim}{numeric vector of length 2, giving limits of response vector to be associated with saturated colors}
  \item{common.lim}{logical, should all plots on a page use the same axis limits}
  \item{...}{any other arguments are passed to \code{imageplot}}
}

\details{
At the time of writing, this function writes plots in PNG format in an arrangement optimized for A4-sized paper.
}

\value{
No value is returned, but one or more files are written to the working directory.
The number of files is determined by the number of columns of \code{RG}.
}

\author{Gordon Smyth}

\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\keyword{hplot}

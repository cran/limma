\name{read.series}
\alias{read.series}
\title{Read series of image files}
\description{
Read in a series of array image analysis output files as data frames.

This function is deprecated in limma 2.9.7 (1 January 2007).
}
\usage{
read.series(slides, path=NULL, suffix="spot", ...)
}
\arguments{
  \item{slides}{character vector giving the names of files to be read in. Any suffix such as ".spot" or ".gpr" which is common to all the files can be omitted.}
  \item{path}{character string giving the directory containing the files. Can be omitted if the files are in the current working directory.}
  \item{suffix}{character string giving a suffix such as "spot" or "gpr" to be added to each file. If \code{NULL} then no suffix is added.}
  \item{\dots}{any other arguments to be passed to \code{read.table}}
}
\details{
This function performs a series of calls to \code{read.table}.
The image analysis output files are assumed to have been edited to remove all pre-heading material.
The files are assumed to contain only column names and data.
In most cases only can use \code{read.maimages} instead.

The data.frames produced by this command will typically be processed further using one of the functions \code{rg.spot}, \code{rg.genepix} or \code{rg.quantarray}.
}
\value{
No value is returned.
However a series of data.frames are created on the current environment with names of the form filename.suffix.
The files names are given by the elements of \code{slides} and the suffix is given by \code{suffix}.
}
\author{Gordon Smyth}
\seealso{
\code{\link[base]{read.table}}.

An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\keyword{file}

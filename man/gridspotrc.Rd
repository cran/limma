\name{gridr}
\alias{gridr}
\alias{gridc}
\alias{spotr}
\alias{spotc}
\title{Row and Column Positions on Microarray}
\description{
Grid and spot row and column positions.
}
\usage{
gridr(layout)
gridc(layout)
spotr(layout)
spotc(layout)
}
\arguments{
  \item{layout}{list with the components \code{ngrid.r}, \code{ngrid.c}, \code{nspot.r} and \code{nspot.c}}
}
\value{
Vector of length \code{prod(unlist(layout))} giving the grid rows (\code{gridr}), grid columns (\code{gridc}), spot rows (\code{spotr}) or spot columns (\code{spotc}).
}
\author{Gordon Smyth}
\keyword{IO}

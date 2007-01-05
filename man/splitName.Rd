\name{splitName}
\alias{splitName}

\title{Split Composite Gene Names}
\description{
Split composite gene names into short names and gene annotation strings.
NOTE: this function is now deprecated in favor of \code{strsplit2}.}

\usage{
splitName(x, split=";", extended=TRUE)
}

\arguments{
\item{x}{character vector}
\item{split}{character to split each element of vector on, see \code{strsplit}}
\item{extended}{logical.  If \code{TRUE}, extended regular expression matching is used, see \code{strsplit}.}
}

\value{
A list containing components
\item{Name}{character vector of the same length as \code{x} contain first splits of each element}
\item{Annotation}{character vector of the same length as \code{x} contain second splits of each element}
}

\details{
Gene names are assumed to comprise a short name or identifier followed by more detailed annotation information.
}

\seealso{
\code{\link[base]{strsplit}}.

An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\author{Gordon Smyth}

\keyword{character}

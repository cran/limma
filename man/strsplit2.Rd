\name{strsplit2}
\alias{strsplit2}

\title{Split Composite Names}
\description{
Split a vector of composite names into a matrix of simple names.}

\usage{
strsplit2(x, split, extended = TRUE, fixed = FALSE, perl = FALSE)
}

\arguments{
\item{x}{character vector}
\item{split}{character to split each element of vector on, see \code{strsplit}}
\item{extended}{logical.  If \code{TRUE}, extended regular expression matching is used, see \code{strsplit}.}
\item{fixed}{logical. If \code{TRUE} match string exactly, otherwise use regular expressions. Has priority over perl and extended.}
\item{perl}{logical. Should perl-compatible regexps be used? Has priority over extended.}
}

\value{
A list containing components
\item{Name}{character vector of the same length as \code{x} contain first splits of each element}
\item{Annotation}{character vector of the same length as \code{x} contain second splits of each element}
}

\details{
This function is the same as \code{\link[base]{strsplit}} except that the output value is a matrix instead of a list.
The first column of the matrix contains the first component from each element of \code{x}, the second column contains the second components etc.
The number of columns is equal to the maximum number of components for any element of \code{x}.
}

\seealso{
\code{\link[base]{strsplit}}.

An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\examples{
x <- c("AA196000;actinin, alpha 3",
"AA464163;acyl-Coenzyme A dehydrogenase, very long chain",
"3E7;W15277;No Annotation")
strsplit2(x,split=";")
}

\author{Gordon Smyth}

\keyword{character}

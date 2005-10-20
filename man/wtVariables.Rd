\name{wtVariables}
\alias{wtVariables}
\title{Read RGList from Image Analysis Output Files}
\description{
  Extracts variable names from weight functions such as \code{wtarea}, \code{wtflags}, \code{wtIgnore.Filter} and user defined functions for fast loading in \code{read.maimages}. 
}	 
\usage{
wtVariables(x, fun)
}
\arguments{
 \item{x}{A character vector of all columns to search against}
 \item{fun}{A weight function to search for any of the possible columns in x}
}
\details{
This is an internally called function by \code{getColClasses} which constructs a colClasses vector for fast loading of only required columns in \code{read.maimages}
}

\value{
Vector of columns required within read.maimages
}
\author{Gordon Smyth}
\seealso{
\code{\link[limma]{getColClasses}}
\code{\link[base]{read.table}}
}
\examples{
searchcols <- c("area","Flags","Ignore Filter")
wtVariables(searchcols, wtarea)
wtVariables(searchcols, wtflags)
wtVariables(searchcols, wtIgnore.Filter)
rm(searchcols) # tidy up
}
\keyword{file}

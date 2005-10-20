\name{getColClasses}
\alias{getColClasses}
\title{Construct colClasses vector for use within read.table}
\description{
  Construct a \code{colClasses} vector, an argument in \code{\link[base]{read.table}} used in \code{\link[limma]{read.maimages}}
}	 
\usage{
getColClasses(cols, ...)
}
\arguments{
 \item{cols}{A character vector of all columns to search against}
 \item{...}{Accepts any character lists, character vectors, or weight functions to match wanted columns against \code{cols}}
}
\details{
This is an internally called function by \code{\link[limma]{read.maimages}} to create a colClasses vector used in \code{read.table} for fast loading of only required columns in \code{read.maimages}
}

\value{
Character. A named vector of classes to be assumed for the columns. 
Possible values are NA (when \code{\link[base]{type.convert}} is used), "NULL" (when the column is skipped)
}
\author{Marcus Davy}
\seealso{\code{\link[limma]{wtVariables}}}
\examples{
allcnames <- c("Block","Column","Row","Name","ID", "F635 Mean","F532 Mean","Flag","Autoflag")
Annotation <- c("Block","Column")
Columns <- list(R="F635 Mean",G="F532 Mean")
Other.columns <- "Autoflag"

getColClasses(allcnames, Annotation, Columns, Other.columns, wtflags)
rm(allcnames, Annotation, Columns, Other.columns) # tidy up
}
\keyword{file}

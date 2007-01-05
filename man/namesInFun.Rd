\name{namesInFun}
\alias{namesInFun}
\title{Find names in a function definition}
\description{
Given a list of names, detect whether these names appear in the definition of a function.
}	 
\usage{
namesInFun(x, fun)
}
\arguments{
 \item{x}{A character vector of names}
 \item{fun}{A function}
}
\details{
This function is used by \code{read.maimages} as part of determining which columns needs to be read from the data file.
It is used to find detect which of the column headings in the data file appear in the definition of the weight function specified by \code{wt.fun}.
}

\value{
Character of vector of names found in the definition of \code{fun}
}
\author{Gordon Smyth}
\seealso{
\code{\link{getColClasses}}
\code{\link{read.maimages}}

An overview of limma read functions is given in \link{03.ReadingData}.
}
\examples{
searchcols <- c("area","Flags","Ignore Filter")
namesInFun(searchcols, wtarea)
namesInFun(searchcols, wtflags)
namesInFun(searchcols, wtIgnore.Filter)
rm(searchcols) # tidy up
}
\keyword{file}

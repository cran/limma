\name{summary}
\alias{summary.RGList}
\alias{summary.MAList}
\alias{summary.MArrayLM}
\title{Summaries of Microarray Data Objects}
\description{
Briefly summarize microarray data objects.
}
\usage{
\method{summary}{RGList}(object, ...)
}
\arguments{
  \item{object}{an object of class \code{RGList}, \code{MAList} or \code{MArrayLM}}
  \item{...}{other arguments are not used}
}
\details{
The data objects are summarized as if they were lists, i.e., brief information about the components is listed.
}
\value{
A table.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{summary}} in the base package.
  
  \link{2.Classes} gives an overview of data classes used in LIMMA.
}
\keyword{methods}

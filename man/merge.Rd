\name{merge}
\alias{merge.RGList}
\alias{merge.MAList}
\title{Merge RGList or MAList Data Objects}
\description{
Merge two microarray data sets represented by RGLists in possibly irregular order.
}
\usage{
\method{merge}{RGList}(x,y,\dots)
}
\arguments{
  \item{x}{\code{\link{RGList-class}} or \code{\link{MAList-class}} object}
  \item{y}{\code{RGList} object, corresponding to the same genes as for \code{x}, possibly in a different order, but with different arrays.}
  \item{\dots}{other arguments are accepted but not used at present}
}
\details{
\code{RGList} and \code{MAList} objects are list objects containing numeric matrices all of the same dimensions.
The RGLists are merged by merging each of the components by row names or, if there are no row names, by IDs in the \code{genes} component.
Unlike when using \code{\link{cbind}}, row names are not required to be in the same order or to be unique.
In the case of repeated row names, the order of the rows with repeated names in preserved.
This means that the first occurrence of each name in \code{x$R} is matched with the first occurrence of the same name in \code{y$R}, the second with the second, and so on.
The final vector of row names is the same as in \code{x}.

Note: if the \code{RGList} objects contain the same number of genes in the same order then the appropriate function to combine them is \code{\link[limma:cbind]{cbind}} rather than \code{merge}.
}
\value{
An merged object of the same class as \code{x} and \code{y} with the same components as \code{x}.
Component matrices have the same rows names as in \code{x} but columns from \code{y} as well as \code{x}.
}
\author{Gordon Smyth}
\seealso{
R base provides a \code{\link[base]{merge}} method for merging data.frames.

An overview of limma commands for reading, subsetting and merging data is given in \link{03.ReadingData}.
}
\examples{
M <- A <- matrix(11:14,4,2)
rownames(M) <- rownames(A) <- c("a","a","b","c")
MA1 <- new("MAList",list(M=M,A=A))

M <- A <- matrix(21:24,4,2)
rownames(M) <- rownames(A) <- c("b","a","a","c")
MA2 <- new("MAList",list(M=M,A=A))

merge(MA1,MA2)
merge(MA2,MA1)
}
\keyword{manip}

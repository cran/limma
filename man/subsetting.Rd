\name{subsetting}
\alias{[.RGList}
\alias{[.MAList}
\alias{[.MArrayLM}
\title{Subset RGList, MAList or MArrayLM Objects}
\description{
Extract a subset of an \code{RGList}, \code{MAList} or \code{MArrayLM} object.
}
\usage{
\method{[}{RGList}(object, i, j, \ldots)
}
\arguments{
  \item{object}{object of class \code{RGList}, \code{MAList} or \code{MArrayLM}}
  \item{i,j}{elements to extract. \code{i} subsets the genes or spots while \code{j} subsets the arrays}
  \item{\ldots}{not used}
}
\details{
\code{i,j} may take any values acceptable for the matrix components of \code{object}.
See the \link[base]{Extract} help entry for more details on subsetting matrices.
}
\value{
An object of the same class as \code{object} holding data from the specified subset of genes and arrays.
}
\author{Gordon Smyth}
\seealso{
  \code{\link[base]{Extract}} in the base package.
  
  \link{03.ReadingData} gives an overview of data input and manipulation functions in LIMMA.
}
\examples{
M <- A <- matrix(11:14,4,2)
rownames(M) <- rownames(A) <- c("a","b","c","d")
colnames(M) <- colnames(A) <- c("A","B")
MA <- new("MAList",list(M=M,A=A))
MA[1:2,]
MA[1:2,2]
MA[,2]
}
\keyword{manip}

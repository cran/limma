\name{bwss.matrix}
\alias{bwss.matrix}
\title{Between and within sums of squares for matrix}
\description{Sums of squares between and within the columns of a matrix. Allows for missing values. This function is called by the \code{\link[limma:anova-method]{anova}} method for \code{MAList} objects.}
\usage{bwss.matrix(x)}
\arguments{
  \item{x}{a numeric matrix.}
}
\value{
  A list with components
  \item{bss}{sums of squares between the column means.}
  \item{wss}{sums of squares within the column means.}
  \item{bdf}{degrees of freedom corresponding to \code{bss}.}
  \item{wdf}{degrees of freedom corresponding to \code{wss}.}
}
\details{This is equivalent to a one-way analysis of variance where the columns of the matrix are the groups.
If \code{x} is a matrix then \code{bwss.matrix(x)} is the same as \code{bwss(x,col(x))} except for speed of execution.}
\author{Gordon Smyth}
\seealso{\code{\link{bwss}}, \code{\link{anova.MAList}}}
\keyword{models}

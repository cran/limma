\name{venn}
\alias{vennCounts}
\alias{vennDiagram}
\title{Venn Diagrams}
\description{
Compute classification counts or plot classification counts in a Venn diagram.
}
\usage{
vennCounts(x, include="both")
vennDiagram(object, include="both", names, mar=rep(1,4), cex=1.5, \dots)
}
\arguments{
  \item{x}{numeric matrix of 0's and 1's indicating significance of a test.
  Usually created by \code{\link{decideTests}}.}
  \item{object}{either a \code{TestResults} matrix or a \code{VennCounts} object produced by \code{vennCounts}.}
  \item{include}{character string specifying whether to counts genes up-regulated, down-regulated or both.
  Choices are \code{"both"}, \code{"up"} or \code{"down"}.}
  \item{names}{optional character vector giving names for the sets or contrasts}
  \item{mar}{numeric vector of length 4 specifying the width of the margins around the plot. This argument is passed to \code{par}.}
  \item{cex}{numerical value giving the amount by which the contrast names should be scaled on the plot relative to the default.plotting text. See \code{par}.}
  \item{\dots}{any other arguments are passed to \code{plot}}
}
\value{
\code{vennCounts} produces a \code{VennCounts} object, which is a numeric matrix with last column \code{"Counts"} giving counts for each possible vector outcome.
\code{vennDiagram} causes a plot to be produced on the current graphical device.
For \code{venDiagram}, the number of columns of \code{object} should be three or fewer.
}
\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\author{Gordon Smyth and James Wettenhall}
\examples{
Y <- matrix(rnorm(100*6),100,6)
Y[1:10,3:4] <- Y[1:10,3:4]+3
Y[1:20,5:6] <- Y[1:20,5:6]+3
design <- cbind(1,c(0,0,1,1,0,0),c(0,0,0,0,1,1))
fit <- eBayes(lmFit(Y,design))
results <- decideTests(fit)
a <- vennCounts(results)
print(a)
vennDiagram(a)
}
\keyword{htest}

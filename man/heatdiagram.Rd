\title{Stemmed Heat Diagram}
\name{heatdiagram}
\alias{heatdiagram}
\alias{heatDiagram}
\description{
Creates a heat diagram showing the co-regulation of genes under one condition with a range of other conditions.
}
\usage{
heatDiagram(results,coef,primary=1,names=NULL,treatments=colnames(coef),limit=NULL,orientation="landscape",cex=1,low="green",high="red",ncolors=123,...)
heatdiagram(stat,coef,primary=1,names=NULL,treatments=colnames(stat),critical.primary=4,critical.other=3,limit=NULL,orientation="landscape",cex=1,low="green",high="red",ncolors=123,...)
}
\arguments{
  \item{results}{\code{TestResults} matrix, containing elements -1, 0 or 1, from \code{\link{classifyTests}}}
  \item{stat}{numeric matrix of test statistics. Rows correspond to genes and columns to treatments or contrasts between treatments.}
  \item{coef}{numeric matrix of the same size as \code{stat}. Holds the coefficients to be displayed in the plot.}
  \item{primary}{number or name of the column to be compared to the others. Genes are included in the diagram according to this column of \code{stat} and are sorted according to this column of \code{coef}. If \code{primary} is a name, then \code{stat} and \code{coef} must have the same column names.}
  \item{names}{optional character vector of gene names}
  \item{treatments}{optional character vector of treatment names}
  \item{critical.primary}{critical value above which the test statistics for the primary column are considered significant and included in the plot}
  \item{critical.other}{critical value above which the other test statistics are considered significant. Should usually be no larger than \code{critical.primary} although larger values are permitted.}
  \item{limit}{optional value for \code{coef} above which values will be plotted in extreme color. Defaults to \code{max(abs(coef))}.}
  \item{orientation}{\code{"portrait"} for upright plot or \code{"landscape"} for plot orientated to be wider than high. \code{"portrait"} is likely to be appropriate for inclusion in printed document while \code{"landscape"} may be appropriate for a presentation on a computer screen.}
  \item{low}{color associated with repressed gene regulation}
  \item{cex}{factor to increase or decrease size of column and row text}
  \item{high}{color associated with induced gene regulation}
  \item{ncolors}{number of distinct colors used for each of up and down regulation}
  \item{...}{any other arguments will be passed to the \code{image} function}
}
\details{
Users are encouraged to use \code{heatDiagram} rather than \code{heatdiagram} as the later function may be removed in future versions of limma.

This function plots an image of gene expression profiles in which rows (or columns for portrait orientation) correspond to treatment conditions and columns (or rows) correspond to genes.
Only genes which are significantly differentially expressed in the primary condition are included.
Genes are sorted by differential expression under the primary condition.

Note: the plot produced by this function is unique to the limma package.
It should not be confused with "heatmaps" often used to display results from cluster analyses.
}
\value{An image is created on the current graphics device.
A dataframe containing the coefficients used in the plot is also invisibly returned.}
\author{Gordon Smyth}
\seealso{\code{\link[graphics]{image}}.}
\examples{
library(sma)
data(MouseArray)
MA <- stat.ma(mouse.data,layout=mouse.setup)
design <- cbind(c(1,1,1,0,0,0),c(0,0,0,1,1,1))
fit <- lm.series(MA$M,design=design)
contrasts.mouse <- cbind(c(1,0),c(0,1),c(-1,1))
colnames(contrasts.mouse) <- c("First3","Second3","Difference")
fit <- contrasts.fit(fit,contrasts=contrasts.mouse)
eb <- ebayes(fit)
heatdiagram(abs(eb$t),fit$coef,primary="Difference")
}
\keyword{hplot}

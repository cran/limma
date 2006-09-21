\title{plotlines}
\name{plotlines}
\alias{plotlines}
\description{
Time course style plot of expression data.
}
\usage{
plotlines(x,first.column.origin=FALSE,xlab="Column",ylab="x",col="black",lwd=1,...)
}
\arguments{
  \item{x}{numeric matrix or object containing expression data.}
  \item{first.column.origin}{logical, should the lines be started from zero?}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{col}{vector of colors for lines}
  \item{lwd}{line width multiplier}
  \item{...}{any other arguments are passed to \code{plot}}
}

\details{
Plots a line for each probe.
}

\value{A plot is created on the current graphics device.}
\author{Gordon Smyth}
\seealso{
An overview of modeling functions and associated plots available in LIMMA is given in \link{06.LinearModels}.
}
\keyword{hplot}

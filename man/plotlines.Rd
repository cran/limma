\title{plotlines}
\name{plotlines}
\alias{plotlines}
\description{
Time course style plot of expression data.
}
\usage{
plotlines(x,first.column.origin=FALSE,xlab="Column",ylab="x",...)
}
\arguments{
  \item{x}{numeric matrix or object containing expression data.}
  \item{first.column.origin}{logical, should the lines be started from zero?}
  \item{xlab}{x-axis label}
  \item{ylab}{y-axis label}
  \item{...}{any other arguments are passed to \code{plot}}
}

\details{
Plots line for each probe.
}

\value{A plot is created on the current graphics device.}
\author{Gordon Smyth}
\seealso{
An overview of diagnostic functions available in LIMMA is given in \link{09.Diagnostics}.
}
\keyword{hplot}

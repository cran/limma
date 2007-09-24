\name{TestResults-class}
\docType{class}
\alias{TestResults-class}
\alias{show,TestResults-method}
\alias{summary.TestResults}
\title{Matrix of Test Results - class}

\description{
A matrix-based class for storing the results of simultanous tests.
\code{TestResults} objects are normally created by \code{\link{classifyTestsF}}, \code{\link{classifyTestsT}} or \code{\link{classifyTestsP}}.
}

\usage{
\S3method{summary}{TestResults}(object, ...)
}

\arguments{
\item{object}{object of class \code{TestResults}}
\item{\dots}{other arguments are not used}
}

\section{Slots/List Components}{
\code{TestResults} objects can be created by \code{new("TestResults",results)} where \code{results} is a matrix.
Objects of this class contain no slots (other than \code{.Data}), although the attributes \code{dim} and \code{dimnames} may be treated as slots.
}

\section{Methods}{
This class inherits directly from class \code{matrix} so any operation appropriate for matrices will work on objects of this class.
\code{\link[methods]{show}} and \code{summary} methods are also implemented.

Functions in LIMMA which operate on \code{TestResults} objects include
\code{\link{heatDiagram}},
\code{\link{vennCounts}},
\code{\link{vennDiagram}},
\code{\link{write.fit}}.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
  \link{08.Tests} gives an overview of multiple testing.
}

\examples{\dontrun{
#  Assume a data object MA and a design matrix
fit <- lmFit(MA, design)
fit <- eBayes(fit)
results <- decideTests(fit)
summary(results)
}}

\keyword{classes}
\keyword{htest}

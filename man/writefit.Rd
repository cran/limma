\name{write.fit}
\alias{write.fit}

\title{Write MArrayLM Object to a File}

\description{
Write a microarray linear model fit a file.
}

\usage{
write.fit(fit, results=NULL, file, digits=2, adjust="none", sep="\t", ...)
}

\arguments{
\item{fit}{object of class \code{MArrayLM} containing the results of a linear model fit}
\item{results}{object of class \code{TestResults}}
\item{file}{character string giving name of file}
\item{digits}{integer indicating precision to be used}
\item{adjust}{character string specifying method to use to adjust the P-values for multiple testing, e.g., "holm" or "fdr". See \code{\link[ctest]{p.adjust}} for the available options. If \code{NULL} or \code{"none"} then the P-values are not adjusted.}
\item{sep}{the field separator string. Values in the output file will be separated by this string.}
\item{\dots}{other arguments are passed to \code{write.table}}
}

\value{
No value is produced but a file is written to the current working directory.

The p-value adjustment is done separately for each column.
This is equivalent to \code{decideTests} with \code{method="separate"}.
}

\details{
This function writes a tab-delimited text file containing for each gene (1) the average log-intensity, (2) the log-ratios, (3) moderated t-statistics, (4) F-statistic if available, (5) classification if available and (6) gene names and annotation.
}

\author{Gordon Smyth}
\seealso{
\code{\link[base]{write}} in the base library.

An overview of linear model functions in limma is given by \link{06.LinearModels}.
}
\keyword{IO}

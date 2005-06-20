\name{08.Tests}
\alias{08.Tests}
\title{Hypothesis Testing for Linear Models}

\description{
LIMMA provides a number of functions for multiple testing across both contrasts and genes.
The starting point is an \code{MArrayLM} object, called \code{fit} say, resulting from fitting a linear model and running \code{eBayes} and, optionally, \code{contrasts.fit}.
See \link{06.LinearModels} or \link{07.SingleChannel} for details.
}

\section{Multiple testing across genes and contrasts}{
The key function is \code{\link{decideTests}}.
This function writes an object of class \code{\link[limma:TestResults]{TestResults}}, which is basically a matrix of \code{-1}, \code{0} or \code{1} elements, of the same dimension as \code{fit$coefficients}, indicating whether each coefficient is significantly different from zero.
A number of different multiple testing strategies are provided.
The function calls other functions \code{\link{classifyTestsF}}, \code{\link{classifyTestsP}} and \code{\link{classifyTestsT}} which implement particular strategies. 
The function \code{\link{FStat}} provides an alternative interface to \code{classifyTestsF} to extract only the overall moderated F-statistic.

A number of other functions are provided to display the results of \code{decideTests}.
The functions \code{\link{heatDiagram}} (or the older version \code{\link{heatdiagram}} displays the results in a heat-map style display.
This allows visual comparison of the results across many different conditions in the linear model.

The functions \code{\link{vennCounts}} and \code{\link{vennDiagram}} provide Venn diagram style summaries of the results.

Summary and \code{show} method exists for objects of class \code{TestResults}.

The results from \code{decideTests} can also be included when the results of a linear model fit are written to a file using \code{\link{write.fit}}.
}

\section{Other Functions}{
Given a set of p-values, the function \code{\link{convest}} can be used to estimate the proportion of true null hypotheses.

When evaluating test procedures with simulated or known results, the utility function \code{\link{auROC}} can be used to compute the area under the Receiver Operating Curve for the test results for a given probe.
}

\author{Gordon Smyth}
\keyword{documentation}

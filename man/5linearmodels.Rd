\name{5.LinearModels}
\alias{5.LinearModels}
\title{Linear Models for Microarrays}

\description{
This page gives an overview of the LIMMA functions available to fit linear models and to interpret the results.
This page covers models for two color arrays in terms of log-ratios or for single-channel arrays in terms of log-intensities.
If you wish to fit models to the individual channel log-intensities from two colour arrays, see \link{6.SingleChannel}.

The core of this package is the fitting of gene-wise linear models to microarray data.
The basic idea is to estimate log-ratios between two or more target RNA samples simultaneously.
See the \emph{\link[limma:../doc/usersguide]{LIMMA User's Guide}} for several case studies.
}

\section{Fitting Models}{

The main function for model fitting is \code{\link{lmFit}}.
This is recommended interface for most users.
\code{lmFit} produces a fitted model object of class \code{\link[limma:marraylm]{MArrayLM}} containing coefficients, standard errors and residual standard errors for each gene.
\code{lmFit} calls one of the following three functions to do the actual computations:

\describe{
\item{ \code{\link{lm.series}} }{
	Straightforward least squares fitting of a linear model for each gene.}

\item{ \code{\link{rlm.series}} }{
	An alternative to \code{lm.series} using robust regression as implemented by the \code{rlm} function in the MASS package.}

\item{ \code{\link{gls.series}} }{
	Generalized least squares taking into account correlations between duplicate spots (i.e., replicate spots on the same array) or related arrays.
	The function \code{\link{duplicateCorrelation}} is used to estimate the inter-duplicate or inter-block correlation before using \code{gls.series}.}
}

All the functions which fit linear models use \code{\link{unwrapdups}} which provides an unified method for handling duplicate spots.
}

\section{Forming the Design Matrix}{

\code{lmFit} has two main arguments, the expression data and the design matrix.
The design matrix is essentially an indicator matrix which specifies which target RNA samples were applied to each channel on each array.
There is considerable freedom in choosing the design matrix - there is always more than one choice which is correct provided it is interpreted correctly.

Design matrices for Affymetrix or single-color arrays can be created using the function \code{\link[stats]{model.matrix}} which is part of the R base package.
The function \code{\link{modelMatrix}} is provided to assist with creation of an appropriate design matrix for two-color microarray experiments.
For direct two-color designs, without a common reference, the design matrix often needs to be created by hand.
}

\section{Making Comparisons of Interest}{

Once a linear model has been fit using an appropriate design matrix, the command \code{\link{makeContrasts}} may be used to form a contrast matrix to make comparisons of interest.
The fit and the contrast matrix are used by \code{\link{contrasts.fit}} to compute fold changes and t-statistics for the contrasts of interest.
This is a way to compute all possible pairwise comparisons between treatments for example in an experiment which compares many treatments to a common reference.
}

\section{Assessing Differential Expression}{

After fitting a linear model, the standard errors are moderated using a simple empirical Bayes model using \code{\link{ebayes}} or \code{\link{eBayes}}.
A moderated t-statistic and a log-odds of differential expression is computed for each contrast for each gene.

\code{\link{ebayes}} and \code{\link{eBayes}} use internal functions \code{\link{fitFDist}}, \code{\link{tmixture.matrix}} and \code{\link{tmixture.vector}}.

The function \code{\link{zscoreT}} is sometimes used for computing z-score equivalents for t-statistics so as to place t-statistics with different degrees of freedom on the same scale.
\code{\link{zscoreGamma}} is used the same way with standard deviations instead of t-statistics.
These functions are for research purposes rather than for routine use.
}

\section{Summarizing Model Fits}{

After the above steps the results may be displayed or further processed using:
\describe{
\item{ \code{\link{toptable}} or \code{\link{topTable}} }{
	Presents a list of the genes most likely to be differentially expressed for a given contrast.}

\item{ \code{\link{classifyTestsF}} }{
	Uses nested F-tests to classify the genes as up, down or even over the contrasts in the linear model with special attention to genes which are significant in more than one contrast.
	\code{\link{classifyTestsT}} and \code{\link{classifyTestsP}} are simpler methods using cutoffs for the t-statistics or p-values individually.}

\item{ \code{\link{FStat}} }{
	Computes an overall moderated F-statistic to test whether all the contrasts are equal to zero.}

\item{ \code{\link{heatdiagram}} or \code{\link{heatDiagram}} }{
	Allows visual comparison of the results across many different conditions in the linear model.
	Not the same as heatdiagrams produced by other packages!
	This function accepts a \code{TestResults} matrix produced by \code{classifyTests}.}

\item{ \code{\link{vennCounts}} }{
	Accepts output from \code{classifyTests} and counts the number of genes in each classification.}

\item{ \code{\link{vennDiagram}} }{
	Accepts output from \code{classifyTests} or \code{vennCounts} and produces a Venn diagram plot.}

\item{ \code{\link{write.fit}} }{
	Writes an \code{MarrayLM} object to a file.}
}
When evaluating test procedures with simulated or known results, the utility function \code{\link{auROC}} can be used to compute the area under the Receiver Operating Curve for the test results for a given probe.
}

\author{Gordon Smyth}
\references{
Smyth, G. K. (2004). Linear models and empirical Bayes methods for assessing differential expression in microarray experiments.
\emph{Statistical Applications in Genetics and Molecular Biology}, \bold{3}, No. 1, Article 3. \url{http://www.bepress.com/sagmb/vol3/iss1/art3}

Smyth, G. K., Michaud, J., and Scott, H. (2003). The use of within-array duplicate spots for assessing differential expression in microarray experiments.
\url{http://www.statsci.org/smyth/pubs/dupcor.pdf}
}
\keyword{documentation}

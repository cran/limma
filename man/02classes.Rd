\name{02.Classes}
\alias{02.Classes}
\title{Classes Defined by this Package}

\description{

This package defines the following data classes.
\describe{
\item{ \code{\link[limma:rglist]{RGList}} }{
	A class used to store raw intensities as they are read in from an image analysis output file,
	usually by \code{read.maimages}.}

\item{ \code{\link[limma:malist]{MAList}} }{
	Intensities converted to M-values and A-values, i.e., to with-spot and whole-spot contrasts on the log-scale.
	Usually created from an \code{RGList} using \code{MA.RG} or \code{normalizeWithinArrays}.
	Objects of this class contain one row for each spot.
	There may be more than one spot and therefore more than one row for each probe.}

\item{ \code{\link[limma:marraylm]{MArrayLM}} }{
	Store the result of fitting gene-wise linear models to the normalized intensities or log-ratios.
	Usually created by \code{lmFit}.
	Objects of this class normally contain only one row for each unique probe.}

\item{ \code{\link[limma:TestResults]{TestResults}} }{
	Store the results of testing a set of contrasts equal to zero for each probe.
	Usually created by \code{decideTests}.
	Objects of this class normally contain one row for each unique probe.}
}

All these data classes obey many analogies with matrices.
In the case of \code{RGList} and \code{MAList}, rows correspond to spots and columns to arrays.
In the case of \code{MarrayLM}, rows correspond to unique probes and the columns to parameters or contrasts.
The functions \code{summary}, \code{dim}, \code{length}, \code{ncol}, \code{nrow}, \code{dimnames}, \code{rownames}, \code{colnames} have methods for these classes.
Objects of any of these classes may be \link[limma:subsetting]{subsetted}.
Multiple data objects may be \link[limma:cbind]{combined} by rows (to add extra probes) or by columns (to add extra arrays).

Furthermore all of these classes may be coerced to actually be of class \code{matrix} using \code{as.matrix}, although this entails loss of information.
Fitted model objects of class \code{MArrayLM} can be coerced to class \code{data.frame} using \code{\link[limma:asdataframe]{as.data.frame}}.

The first three classes belong to the virtual class \code{\link[limma:LargeDataObject]{LargeDataObject}}.
A \code{show} method is defined for \code{LargeDataOject}s which uses the utility function \code{\link{printHead}}.
}

\author{Gordon Smyth}
\keyword{documentation}

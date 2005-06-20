\name{exprSet2-class}
\docType{class}
\alias{exprSet2-class}
\title{Expression Set - class}

\description{
A class for storing intensity values from microarray experiments.
This class is similar to \code{\link[Biobase:exprSet-class]{exprSet}}.

This class is not currently used in the limma package.
It is an experimental class to investigate the possibility of a unified data class structure for single-channel and log-ratio analyses of spotted microarray data.
}

\section{Slots}{
\describe{
	\item{expressions}{\code{matrix} containing intensity data on log-2 scale.
	For two-color arrays, odd columns will usually correspond to channel 1 (green) and even columns to channel 2 (red) for different arrays.}
	\item{weights}{\code{matrix} containing non-negative quality weights}
	\item{targets}{\code{data.frame} containing factors corresponding to the columns of \code{expressions}}
	\item{probes}{\code{data.frame} containing gene IDs or annotation information.
	Should have same number of rows as \code{expression}}
	\item{printer}{\code{list} containing information about the printing process}
	\item{notes}{\code{character}}
}
}

\section{Methods}{
\code{exprSet2} objects inherit a \code{\link[methods]{show}} method from the virtual class \code{\link[limma:LargeDataObject]{LargeDataObject}}, which means that \code{exprSet2} objects will print in a compact way.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
  
  \code{\link[Biobase:exprSet-class]{exprSet}} is the corresponding class in the Biobase package.
}

\keyword{classes}
\keyword{data}

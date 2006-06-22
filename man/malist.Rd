\name{MAList-class}
\docType{class}
\alias{MAList-class}
\title{M-value, A-value Expression List - class}

\description{
A simple list-based class for storing M-values and A-values for a batch of spotted microarrays.
\code{MAList} objects are usually created during normalization by the functions \code{\link{normalizeWithinArrays}} or \code{\link{MA.RG}}.
}

\section{Slots/List Components}{
\code{MAList} objects can be created by \code{new("MAList",MA)} where \code{MA} is a list.
This class contains no slots (other than \code{.Data}), but objects should contain the following components:
\tabular{ll}{
  \code{M}:\tab numeric matrix containing the M-values (log-2 expression ratios).  Rows correspond to spots and columns to arrays.\cr
  \code{A}:\tab numeric matrix containing the A-values (average log-2 expression values).
}
Optional components include:
\tabular{ll}{
  \code{weights}:\tab numeric matrix of same dimensions as \code{M} containing relative spot quality weights.  Elements should be non-negative.\cr
  \code{other}:\tab list containing numeric matrices of other spot-specific information. All matrices must have the same dimensions as \code{M}.\cr
  \code{genes}:\tab data.frame containing probe information. Should have one row for each spot. May have any number of columns.\cr
  \code{targets}:\tab data.frame containing information on the target RNA samples.  Rows correspond to arrays.  May have any number of columns. Usually includes columns \code{Cy3} and \code{Cy5} specifying which RNA was hybridized to each array.\cr
  \code{printer}:\tab list containing information on the process used to print the spots on the arrays.  See \link[limma:PrintLayout]{PrintLayout}.
}
Valid \code{MAList} objects may contain other optional components, but all probe or array information should be contained in the above components.
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
In addition, \code{MAList} objects can be \link[limma:subsetting]{subsetted} and \link[limma:cbind]{combined}.
\code{RGList} objects will return dimensions and hence functions such as \code{\link[limma:dim]{dim}}, \code{\link[base:nrow]{nrow}} and \code{\link[base:nrow]{ncol}} are defined. 
\code{MALists} also inherit a \code{\link[methods]{show}} method from the virtual class \code{\link[limma:LargeDataObject]{LargeDataObject}}, which means that \code{RGLists} will print in a compact way.

Other functions in LIMMA which operate on \code{MAList} objects include
\code{\link{normalizeWithinArrays}},
\code{\link{normalizeBetweenArrays}},
\code{\link{normalizeForPrintorder}},
\code{\link{plotMA}}
and \code{\link{plotPrintTipLoess}}.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
  
  \code{\link[marrayClasses]{marrayNorm-class}} is the corresponding class in the marrayClasses package.
}

\keyword{classes}
\keyword{data}

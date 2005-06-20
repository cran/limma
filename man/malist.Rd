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
This class contains no slots (other than \code{.Data}), but objects should contain the following list components:
\tabular{ll}{
  \code{M}:\tab numeric matrix containing the M-values or log-2 expression ratios.  Rows correspond to spots and columns to arrays.\cr
  \code{A}:\tab numeric matrix containing the A-values or average log-2 expression values
}
\tabular{ll}{
  \code{weights}:\tab numeric matric containing relative spot quality weights.  Should be non-negative.\cr
  \code{printer}:\tab list containing information on the process used to print the spots on the arrays.  See \link[limma:PrintLayout]{PrintLayout}.\cr
  \code{genes}:\tab data.frame containing information on the genes spotted on the arrays.  Should include a character column \code{Name} containing names for the genes or controls.\cr
  \code{targets}:\tab data.frame containing information on the target RNA samples.  Should include factor or character columns \code{Cy3} and \code{Cy5} specifying which RNA was hybridized to each array.
}
All of the matrices should have the same dimensions.
The row dimension of \code{targets} should match the column dimension of the matrices.
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

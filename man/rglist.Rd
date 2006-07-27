\name{RGList-class}
\docType{class}
\alias{RGList-class}
\alias{coerce,RGList,exprSet2-method}
\title{Red, Green Intensity List - class}

\description{
A simple list-based class for storing red and green channel foreground and background intensities for a batch of spotted microarrays.
\code{RGList} objects are normally created by \code{\link{read.maimages}}.
}

\section{Slots/List Components}{
\code{RGList} objects can be created by \code{new("RGList",RG)} where \code{RG} is a list.
Objects of this class contains no slots (other than \code{.Data}), but objects should contain the following list components:
\tabular{ll}{
  \code{R}:\tab numeric matrix containing the red (cy5) foreground intensities.  Rows correspond to spots and columns to arrays.\cr
  \code{G}:\tab numeric matrix containing the green (cy3) foreground intensities. Rows correspond to spots and columns to arrays.
}
Optional components include
\tabular{ll}{
  \code{Rb}:\tab numeric matrix containing the red (cy5) background intensities\cr
  \code{Gb}:\tab numeric matrix containing the green (cy3) background intensities\cr
  \code{weights}:\tab numeric matrix of same dimension as \code{R} containing relative spot quality weights.  Elements should be non-negative.\cr
  \code{other}:\tab list containing numeric matrices of other spot-specific information.\cr
  \code{genes}:\tab data.frame containing probe information. Should have one row for each spot. May have any number of columns.\cr
  \code{targets}:\tab data.frame containing information on the target RNA samples.  Rows correspond to arrays.  May have any number of columns. Usually includes columns \code{Cy3} and \code{Cy5} specifying which RNA was hybridized to each array.\cr
  \code{printer}:\tab list containing information on the process used to print the spots on the arrays.  See \link[limma:PrintLayout]{PrintLayout}.
}
Valid \code{RGList} objects may contain other optional components, but all probe or array information should be contained in the above components.
}

\section{Methods}{
This class inherits directly from class \code{list} so any operation appropriate for lists will work on objects of this class.
In addition, \code{RGList} objects can be \link[limma:subsetting]{subsetted}, \link[limma:cbind]{combined} and \link[limma:merge]{merged}.
\code{RGList} objects will return dimensions and hence functions such as \code{\link[limma:dim]{dim}}, \code{\link[base:nrow]{nrow}} and \code{\link[base:nrow]{ncol}} are defined. 
\code{RGLists} also inherit a \code{\link[methods]{show}} method from the virtual class \code{\link[limma:LargeDataObject]{LargeDataObject}}, which means that \code{RGLists} will print in a compact way.

\code{RGList} objects can be converted to \code{exprSet2} objects by \code{as(RG,"exprSet2")}.

Other functions in LIMMA which operate on \code{RGList} objects include
\code{\link{normalizeBetweenArrays}},
\code{\link{normalizeForPrintorder}},
\code{\link{normalizeWithinArrays}}.
}

\author{Gordon Smyth}

\seealso{
  \link{02.Classes} gives an overview of all the classes defined by this package.
  
  \code{\link[marrayClasses]{marrayRaw-class}} is the corresponding class in the marrayClasses package.
}

\keyword{classes}
\keyword{data}

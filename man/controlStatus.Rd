\name{controlStatus}
\alias{controlStatus}
\title{Set Status of each Spot from List of Spot Types}
\description{
Determine the type (or status) of each spot in the gene list.
}
\usage{
controlStatus(types, genes, spottypecol="SpotType", regexpcol, verbose=TRUE)
}
\arguments{
  \item{types}{dataframe containing spot type specifiers, usually input using \code{readSpotTypes}}
  \item{genes}{dataframe containing the microarray gene list, or an \code{RGList}, \code{MAList} or \code{MArrayList} containing \code{genes}}
  \item{spottypecol}{integer or name specifying column of \code{types} containing spot type names}
  \item{regexpcol}{vector of integers or column names specifying columns of types containing regular expressions.
  Defaults to any column names in common between \code{types} and \code{genes}.}
  \item{verbose}{logical, if \code{TRUE} then progess on pattern matching is reported to the standard output channel}
}
\details{
This function constructs a vector of status codes by searching for patterns in the gene list.
The data frame \code{genes} contains gene IDs and should have as many rows as there are spots on the microarrays.
Such a data frame is often read using \code{\link{readGAL}}.
The data frame \code{types} has as many rows as you want to distinguish types of spots in the gene list.
This data frame should contain a column or columns, the \code{regexpcol} columns, which have the same names as columns in \code{genes} and which contain patterns to match in the gene list.
Another column, the \code{spottypecol}, contains the names of the spot types.
Any other columns are assumed to contain plotting parameters, such as colors or symbols, to be associated with the spot types.

The patterns in the \code{regexpcol} columns are simplified regular expressions.
For example, \code{AA*} means any string starting with \code{AA}, \code{*AA} means any code ending with \code{AA}, \code{AA} means exactly these two letters, \code{*AA*} means any string containing \code{AA}, \code{AA.} means \code{AA} followed by exactly one other character and \code{AA\.} means exactly \code{AA} followed by a period and no other characters.
Any other regular expressions are allowed but the codes \code{^} for beginning of string and \code{$} for end of string should not be included.

Note that the patterns are matched sequentially from first to last, so more general patterns should be included first.
For example, it is often a good idea to include a default spot-type as the first line in \code{types} with pattern \code{*} for all \code{regexpcol} columns and default plotting parameters.
}
\value{
Character vector specifying the type (or status) of each spot on the array.
Attributes contain plotting parameters associated with each spot type.
}
\author{Gordon Smyth}
\seealso{
An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\examples{
genes <- data.frame(ID=c("Control","Control","Control","Control","AA1","AA2","AA3","AA4"),
Name=c("Ratio 1","Ratio 2","House keeping 1","House keeping 2","Gene 1","Gene 2","Gene 3","Gene 4"))
types <- data.frame(SpotType=c("Gene","Ratio","Housekeeping"),ID=c("*","Control","Control"),Name=c("*","Ratio*","House keeping*"),col=c("black","red","blue"))
status <- controlStatus(types,genes)
}
\keyword{IO}

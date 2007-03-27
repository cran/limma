\name{read.maimages}
\alias{read.maimages}
\alias{read.imagene}
\title{Read RGList from Image Analysis Output Files}
\description{
Reads an RGList from a series of two-color microarray image analysis output files
}
\usage{
read.maimages(files=NULL,source="generic",path=NULL,ext=NULL,names=NULL,columns=NULL,other.columns=NULL,annotation=NULL,wt.fun=NULL,verbose=TRUE,sep="\t",quote=NULL,\dots)
read.imagene(files,path=NULL,ext=NULL,names=NULL,columns=NULL,other.columns=NULL,wt.fun=NULL,verbose=TRUE,sep="\t",quote="\"",\dots)
}
\arguments{
  \item{files}{character vector giving the names of the files containing image analysis output or, for Imagene data, a character matrix of names of files.
  If omitted, then all files with extension \code{ext} in the specified directory will be read in alphabetical order.}
  \item{source}{character string specifying the image analysis program which produced the output files.  Choices are \code{"generic"}, \code{"agilent"}, \code{"arrayvision"}, \code{"bluefuse"}, \code{"genepix"}, \code{"genepix.custom"}, \code{"genepix.median"}, \code{"imagene"}, \code{"quantarray"}, \code{"scanarrayexpress"}, \code{"smd.old"}, \code{"smd"}, \code{"spot"} or \code{"spot.close.open"}.}
  \item{path}{character string giving the directory containing the files.
  The default is the current working directory.}
  \item{ext}{character string giving optional extension to be added to each file name}
  \item{names}{character vector of names to be associated with each array as column name.
  Defaults to \code{removeExt(files)}.}
  \item{columns}{list with fields \code{R}, \code{G}, \code{Rb} and \code{Gb} giving the column names to be used for red and green foreground and background or, in the case of Imagene data, a list with fields \code{f} and \code{b}.  This argument is optional if \code{source} is specified, otherwise it is required.}
  \item{other.columns}{character vector of names of other columns to be read containing spot-specific information}
  \item{annotation}{character vector of names of columns containing annotation information about the probes}
  \item{wt.fun}{function to calculate spot quality weights}
  \item{verbose}{logical, \code{TRUE} to report each time a file is read}
  \item{sep}{the field separator character}
  \item{quote}{character string of characters to be treated as quote marks}
  \item{\dots}{any other arguments are passed to \code{read.table}}
}
\details{
This is the main data input function for the LIMMA package for two-color microarray data.
It extracts the foreground and background intensities from a series of files, produced by an image analysis program, and assembles them into the components of one list.
The image analysis programs Agilent Feature Extraction, ArrayVision, BlueFuse, GenePix, ImaGene, QuantArray (Version 3 or later), Stanford Microarray Database (SMD) and SPOT are supported explicitly.
Data from some other image analysis programs can be read if the appropriate column names containing the foreground and background intensities are specified using the \code{columns} argument.
(This will work if the column names are unique and if there are no rows in the file after the last line of data.
Header lines are ok.)

SMD data should consist of raw data files from the database, in tab-delimited text form.
There are two possible sets of column names depending on whether the data was entered into the database before or after September 2003.
\code{source="smd.old"} indicates that column headings in use prior to September 2003 should be used.
In the case of GenePix, two possible foreground estimators are supported: \code{source="genepix"} uses the mean foreground estimates while \code{source="genepix.median"} uses median foreground estimates.
GenePix 6.0 and later also supplies some custom background options, notably morphological background.
If the GPR files have been written using a custom background, you may read it using \code{source="genepix.custom"}. 
In the case of SPOT, two possible background estimators are supported:
if \code{source="spot.close.open"} then background intensities are estimated from \code{morph.close.open} rather than \code{morph}.

Spot quality weights may be extracted from the image analysis files using a weight function wt.fun.
\code{wt.fun} may be any user-supplied function which accepts a data.frame argument and returns a vector of non-negative weights.
The columns of the data.frame are as in the image analysis output files.
There is one restriction, which is that the column names should be refered to in full form in the weight function, i.e., do not rely on name expansion for partial matches when refering to the names of the columns.
See \code{\link{QualityWeights}} for suggested weight functions.

For Imagene image data the argument \code{files} should be a matrix with two columns.
The first column should contain the names of the files containing green channel (cy3) data and the second column should contain names of files containing red channel (cy5) data.
If \code{source="imagene"} and \code{files} is a vector of even length instead of a matrix, then each consecutive pair of file names is assumed to correspond to the same array.
The function \code{read.imagene} is called by \code{read.maimages} when \code{source="imagene"}.
It does not need to be called directly by users.

The argument \code{other.columns} allows arbitrary columns of the image analysis output files to be preserved in the data object.
These become matrices in the component \code{other} component.
For ImaGene data, the other column headings with be prefixed with \code{"R "} or \code{"G "} as appropriate.
}

\section{Warnings}{
All image analysis files being read are assumed to contain data for the same genelist in the same order.
No checking is done to confirm that this is true.
Probe annotation information is read from the first file only.
}

\value{
An \code{\link[limma:rglist]{RGList}} object containing the components
  \item{R}{matrix containing the red channel foreground intensities for each spot for each array.}
  \item{Rb}{matrix containing the red channel background intensities for each spot for each array.}
  \item{G}{matrix containing the green channel foreground intensities for each spot for each array.}
  \item{Gb}{matrix containing the green channel background intensities for each spot for each array.}
  \item{weights}{spot quality weights, if \code{wt.fun} is given}
  \item{other}{list containing matrices corresponding to \code{other.columns} if given}
  \item{genes}{data frame containing annotation information about the probes, for example gene names and IDs and spatial positions on the array, currently set only if \code{source} is \code{"agilent"}, \code{"genepix"} or \code{source="imagene"} or if the \code{annotation} argument is set}
  \item{targets}{data frame with column \code{FileName} giving the names of the files read}
  \item{source}{character string giving the image analysis program name}
  \item{printer}{list of class \code{\link[limma:printlayout]{PrintLayout}}, currently set only if \code{source="imagene"}}
}
\author{Gordon Smyth, with speed improvements by Marcus Davy}
\references{
Web pages for the image analysis software packages mentioned here are listed at \url{http://www.statsci.org/micrarra/image.html}
}
\seealso{
\code{read.maimages} is based on \code{\link[base]{read.table}} in the base package.
\code{\link[marray]{read.marrayRaw}} is the corresponding function in the marray package.

An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}
\examples{
#  Read all .gpr files from current working directory
#  and give weight 0.1 to spots with negative flags

\dontrun{files <- dir(pattern="*\\\\.gpr$")
RG <- read.maimages(files,"genepix",wt.fun=wtflags(0.1))}

#  Read all .spot files from current working director and down-weight
#  spots smaller or larger than 150 pixels

\dontrun{files <- dir(pattern="*\\\\.spot$")
RG <- read.maimages(files,"spot",wt.fun=wtarea(150))}
}
\keyword{file}

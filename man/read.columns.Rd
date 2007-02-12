\name{read.columns}
\alias{read.columns}
\title{Read specified columns from a file}
\description{
Reads specified columns from a file in table format and creates a data frame from it, with cases corresponding to lines and variables to fields in the file. 
}
\usage{
read.columns(file,required.col=NULL,text.to.search="",sep="\t",quote="\"",skip=0,fill=TRUE,blank.lines.skip=TRUE,comment.char="",allowEscapes=FALSE,...)
}
\arguments{
  \item{file}{the name of the file which the data are to be read from.}
  \item{required.col}{character vector of names of the required columns}
  \item{text.to.search}{character string. If any column names can be found in this string, those columns will also be read.}
  \item{sep}{the field separator character}
  \item{quote}{character string of characters to be treated as quote marks}
  \item{skip}{the number of lines of the data file to skip before beginning to read data.} 
  \item{fill}{logical: if \code{TRUE} then in case the rows have unequal length, blank fields are implicitly added.} 
  \item{blank.lines.skip}{logical: if \code{TRUE} blank lines in the input are ignored.}
  \item{comment.char}{character: a character vector of length one containing a single character or an empty string.}
  \item{allowEscapes}{logical. Should C-style escapes such as \code{\n} be processed or read verbatim (the default)?} 
  \item{\dots}{any other arguments are passed to \code{read.table}}
}

\details{
This function is an interface to \code{read.table} in the base package.
It uses \code{required.col} and \code{text.to.search} to set up the \code{colClasses} argument of \code{read.table}.

This function is used by \code{\link{read.maimages}}.
}

\value{
A data frame (data.frame) containing a representation of the data in the file. 
}

\author{Gordon Smyth}

\seealso{
\code{\link{read.maimages}}, \code{\link[base]{read.table}}.

An overview of LIMMA functions for reading data is given in \link{03.ReadingData}.
}

\keyword{file}

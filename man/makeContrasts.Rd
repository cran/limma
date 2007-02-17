\name{makeContrasts}
\alias{makeContrasts}
\title{Construct Matrix of Custom Contrasts}
\description{
Construct the contrast matrix corresponding to specified contrasts of a set of parameters.
}
\usage{
makeContrasts(\dots, contrasts=NULL, levels) 
}
\arguments{
  \item{\dots}{expressions, or character strings which can be parsed to expressions, specifying contrasts}
  \item{contrasts}{character vector specifying contrasts}
  \item{levels}{character vector or factor giving the names of the parameters of which contrasts are desired, or a design matrix or other object with the parameter names as column names.}
}

\value{
Matrix which columns corresponding to contrasts.
}

\details{
This function expresses contrasts between a set of parameters as a numeric matrix.
The parameters are usually the coefficients from a linear model fit, so the matrix specifies which comparisons between the coefficients are to be extracted from the fit.
The output from this function is usually used as input to \code{\link{contrasts.fit}}.
The contrasts can be specified either as expressions using \code{\dots} or as a character vector through \code{contrasts}.
(Trying to specify contrasts both ways will cause an error.)

The parameter names must be syntactically valid variable names in R and so, for example, must begin with a letter rather than a numeral.
See \code{\link{make.names}} for a complete specification of what is a valid name.
}

\seealso{
An overview of linear model functions in limma is given by the help page \link{06.LinearModels}.
}

\author{Gordon Smyth}

\examples{
makeContrasts(B-A,C-B,C-A,levels=c("A","B","C"))
makeContrasts(contrasts="A-(B+C)/2",levels=c("A","B","C"))
x <- c("B-A","C-B","C-A")
makeContrasts(contrasts=x,levels=c("A","B","C"))
}

\keyword{regression}

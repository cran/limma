\name{makeContrasts}
\alias{makeContrasts}
\title{Construct Matrix of Custom Contrasts}
\description{
Construct the contrasts matrix corresponding to specified contrasts of a set of given parameters.
}
\usage{
makeContrasts(\dots, levels) 
}
\arguments{
  \item{\dots}{expressions, or character strings which can be parsed to expressions, specifying contrasts}
  \item{levels}{character vector giving the names of the parameters to be contrasts, or a factor or design matrix from which the names can be extracted.}
}

\value{
Matrix which columns corresponding to contrasts.
}

\seealso{
An overview of linear model functions in limma is given by \link{06.LinearModels}.
}

\author{Gordon Smyth}

\examples{
makeContrasts(B-A,C-B,C-A,levels=c("A","B","C"))
makeContrasts("A","B","B-A",levels=c("A","B"))
}

\keyword{regression}

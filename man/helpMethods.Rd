\name{helpMethods}
\alias{helpMethods}
\title{Prompt for Method Help Topics}
\description{
For any S4 generic function, find all methods defined in currently loaded packages.
Prompt the user to choose one of these to display the help document.
}
\usage{
helpMethods(genericFunction)
}
\arguments{
  \item{genericFunction}{a generic function or a character string giving the name of a generic function}
}
\author{Gordon Smyth}
\seealso{
\code{\link[methods]{showMethods}}
}
\examples{
\dontrun{helpMethods(show)}
}
\keyword{methods}

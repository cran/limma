\name{limmaUsersGuide}
\alias{limmaUsersGuide}

\title{View Limma User's Guide}
\description{Finds the location of the Limma User's Guide and optionally opens it.}

\usage{
limmaUsersGuide(view=TRUE)
}

\arguments{
  \item{view}{logical, should the document be opened using the default PDF document reader?}
}

\value{Character string giving the file location.}

\details{
The function \code{vignette("limma")} will find the short limma Vignette which describes how to obtain the Limma User's Guide.
The User's Guide is not itself a true vignette because it is not automatically generated using \code{\link[utils]{Sweave}} during the package build process.
This means that it cannot be found using \code{vignette}, hence the need for this special function.

If the operating system is other than Windows, then the PDF viewer used is that given by \code{Sys.getenv("R_PDFVIEWER")}.
The PDF viewer can be changed using \code{Sys.putenv(R_PDFVIEWER=)}.

This function is used by drop-down Vignettes menu when the Rgui interface for Windows is used.
}

\seealso{
\code{\link[utils]{vignette}}, \code{\link[Biobase]{openPDF}}, \code{\link[Biobase]{openVignette}}, \code{\link[base]{Sys.getenv}}, \code{\link[base]{Sys.putenv}}
}

\author{Gordon Smyth}

\examples{
limmaUsersGuide(view=FALSE)
}

\keyword{documentation}

\name{uniquegenelist}
\alias{uniquegenelist}
\title{Eliminate Duplicate Names from the Gene List}
\description{
Eliminate duplicate names from the gene list. The new list is shorter than the full list by a factor of \code{ndups}.
}
\usage{
uniquegenelist(genelist,ndups=2,spacing=1)
}
\arguments{
  \item{genelist}{vector of gene names}
  \item{ndups}{number of duplicate spots. The number of rows of \code{genelist} must be divisible by \code{ndups}.}
  \item{spacing}{the spacing between duplicate names in \code{genelist}}
}
\value{
A vector of length \code{length(genelist)/ndups} containing each gene name once only.
}
\author{Gordon Smyth}
\seealso{
\code{\link{unwrapdups}}
}
\examples{
genelist <- c("A","A","B","B","C","C","D","D")
uniquegenelist(genelist,ndups=2)
genelist <- c("A","B","A","B","C","D","C","D")
uniquegenelist(genelist,ndups=2,spacing=2)
}
\keyword{array}

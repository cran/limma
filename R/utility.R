#  UTILITY FUNCTIONS

.matvec <- function(M,v) {
#	Multiply the columns of matrix by the elements of a vector,
#	i.e., compute M %*% diag(v)
#	Gordon Smyth
#	5 July 1999
#
	v <- as.vector(v)
	M <- as.matrix(M)
	if(length(v)!=dim(M)[2]) stop("Dimensions do not match")
	t(v * t(M))
}

.vecmat <- function(v,M) {
#	Multiply the rows of matrix by the elements of a vector,
#	i.e., compute diag(v) %*% M
#	Gordon Smyth
#	5 July 1999
#
	v <- as.vector(v)
	M <- as.matrix(M)
	if(length(v)!=dim(M)[1]) stop("Dimensions do not match")
	v * M
}

isNumeric <- function(x) {
#	Test for numeric argument or data.frame with numeric columns
#	Gordon Smyth
#	12 April 2003

	is.numeric(x) || (is.data.frame(x) && length(x)>0 && all(unlist(lapply(x,is.numeric))))
}

blockDiag <- function(...)
#	Block diagonal matrix
#	Gordon Smyth
#	8 Feb 2004
{
	e <- list(...)
	d <- matrix(unlist(lapply(e,dim)),ncol=2,byrow=TRUE)
	if(nrow(d) != length(e)) stop("all arguments must be matrices")
	dsum <- apply(d,2,sum)
	z <- array(0,dsum)
	dimnames(z) <- list(character(dsum[1]),character(dsum[2]))
	coord <- c(0,0)
	for (i in 1:length(e)) {
		coord1 <- coord[1]+1:d[i,1]
		coord2 <- coord[2]+1:d[i,2]
		z[coord1,coord2] <- e[[i]]
		rownames(z)[coord1] <- rownames(e[[i]],do.NULL=FALSE,prefix="")
		colnames(z)[coord2] <- colnames(e[[i]],do.NULL=FALSE,prefix="")
		coord <- coord+d[i,]
	}
	z
}

helpMethods <- function(genericFunction) {
#	Prompt user for help topics on methods for generic function
#	Gordon Smyth
#	21 April 2003.  Last revised 28 Oct 2003.

	objectclass <- class(genericFunction)
 	if(objectclass != "standardGeneric") {
		if(objectclass == "character" && isGeneric(genericFunction))
			genericFunction <- getGeneric(genericFunction)
		else {
			cat("Not a generic function\n")
			return(invisible())
		}
	}
	functionname <- genericFunction@generic
	methodnames <- names(getMethods(genericFunction)@methods)
	nmethods <- length(methodnames)
	if(nmethods == 0) {
		cat("No available methods\n")
		return(invisible())
	}
	aliasnames <- paste(functionname,",",methodnames,"-method",sep="")
	for (i in 1:nmethods) cat(i,": ",aliasnames[i],"\n",sep="")
	cat("Type number to choose help topic: ")
	n <- as.integer(readline())
	if(n > 0 && n <= nmethods)
		eval(parse(text=paste("help(\"",aliasnames[n],"\")",sep="")))
	else {
	 	cat("No topic chosen\n")
	 	return(invisible())
	}
}

limmaUsersGuide <- function(view=TRUE)
#	Find and optionally view limma User's Guide
#	Gordon Smyth
#	25 Oct 2004.
{
	f <- system.file("doc","usersguide.pdf",package="limma")
	if(view) {
		if(.Platform$OS.type == "windows") 
			shell.exec(f)
		else
			system(paste(Sys.getenv("R_PDFVIEWER"),f,"&"))
	}
	return(f)
}

changeLog <- function(n=20)
#	Write first n lines of limma changelog
#	Gordon Smyth
#	20 Sep 2004.  Last modified 30 Oct 2004.
{
	writeLines(readLines(system.file("doc","changelog.txt",package="limma"),n=n))
}


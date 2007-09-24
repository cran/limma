#	CLASSES.R

setClass("RGList",
#  Class to hold initial read-in data
representation("list")
)

setClass("MAList",
#  Class to hold normalized, rotated data
representation("list")
)

setClass("MArrayLM",
#  Linear model fit
representation("list")
)

setClass("exprSet2",representation(
	expressions="matrix",
	weights="matrix",
	targets="data.frame",
	probes="data.frame",
	printer="list",
	notes="character"
))

setAs("RGList", "exprSet2", function(from, to) {
	y <- new(to)
	if(length(from$G)) y@expressions <- cbind(from$G,from$R)
	if(length(from$weights)) y@weights <- cbind(from$weights,from$weights)
	if(length(from$genes)) y@probes <- from$genes
	if(length(from$printer)) y@printer <- unclass(from$printer)
	y
})

printHead <- function(x)
#  Print leading 5 elements or rows of atomic object
#  Gordon Smyth
#  May 2003.  Last modified 7 April 2004.
{
	if(is.atomic(x)) {
		d <- dim(x)
		if(length(d)<2) which <- "OneD"
		if(length(d)==2) which <- "TwoD"
		if(length(d)>2) which <- "Array"
	} else {
		if(inherits(x,"data.frame")) {
			d <- dim(x)
			which <- "TwoD"
		} else
			which <- "Recursive"
	}
	switch(which,
	OneD={
		n <- length(x)
		if(n > 20) {
			print(x[1:5])
			cat(n-5,"more elements ...\n")
		} else
			print(x)
	},
	TwoD={
		n <- d[1]
		if(n > 10) {
			print(x[1:5,])
			cat(n-5,"more rows ...\n")
		} else
			print(x)
	},
	Array={
		n <- d[1]
		if(n > 10) {
			dn <- dimnames(x)
			dim(x) <- c(d[1],prod(d[-1]))
			x <- x[1:5,]
			dim(x) <- c(5,d[-1])
			if(!is.null(dn[[1]])) dn[[1]] <- dn[[1]][1:5]
			dimnames(x) <- dn
			print(x)
			cat(n-5,"more rows ...\n")
		} else
			print(x)
	},
	Recursive=print(x)
	)
}

setClass("LargeDataObject")
setIs("RGList","LargeDataObject")
setIs("MAList","LargeDataObject")
setIs("MArrayLM","LargeDataObject")
setIs("exprSet2","LargeDataObject")

setMethod("show","LargeDataObject",
#  Print and show method large data objects
#  Gordon Smyth
#  May 2003
function(object) {
	cat("An object of class \"",class(object),"\"\n",sep="")
	for (what in names(object)) {
		x <- object[[what]]
		cat("$",what,"\n",sep="")
		printHead(x)
		cat("\n")
	}
	for (what in setdiff(slotNames(object),".Data")) {
		x <- slot(object,what)
		if(length(x) > 0) {
			cat("@",what,"\n",sep="")
			printHead(x)
			cat("\n")
		}
	}
})

dim.RGList <- function(x) if(is.null(x$R)) c(0,0) else dim(as.matrix(x$R))
dim.MAList <- function(x) if(is.null(x$M)) c(0,0) else dim(as.matrix(x$M))
dim.MArrayLM <- function(x) if(is.null(x$coefficients)) c(0,0) else dim(as.matrix(x$coefficients))
length.RGList <- length.MAList <- length.MArrayLM <- function(x) prod(dim(x))

dimnames.RGList <- function(x) dimnames(x$R)
dimnames.MAList <- function(x) dimnames(x$M)
dimnames.MArrayLM <- function(x) dimnames(x$coefficients)
.setdimnames <- function(x, value)
#  Dimension names for RGList-like objects
#  Gordon Smyth
#  17 Dec 2005.
{
	exprmatrices <- c("R","G","Rb","Gb","M","A","weights")
	for (a in exprmatrices) if(!is.null(x[[a]])) dimnames(x[[a]]) <- value
	for(a in names(x$other)) dimnames(object$other[[a]]) <- value
	if(!is.null(x$targets)) row.names(x$targets) <- value[[2]]
	if(!is.null(x$design)) rownames(x$design) <- value[[2]]
	x
}
#assign("dimnames<-.RGList",.setdimnames)
#assign("dimnames<-.MAList",.setdimnames)
"dimnames<-.RGList" <- .setdimnames
"dimnames<-.MAList" <- .setdimnames

summary.MArrayLM <- summary.MAList <- summary.RGList <- function(object,...) summary(unclass(object))

as.MAList <- function(object) {
#	Convert marrayNorm object to MAList
#	Gordon Smyth
#	20 Sep 2003.  Last modified 20 Dec 2003.

	MA <- new("MAList")
	ifposlen <- function(x) if(length(x)) return(x) else return(NULL)
	MA$A <- ifposlen(object@maA)
	MA$M <- ifposlen(object@maM)
	MA$weights <- ifposlen(object@maW)
	MA$printer$ngrid.r <- ifposlen(object@maLayout@maNgr)
	MA$printer$ngrid.c <- ifposlen(object@maLayout@maNgc)
	MA$printer$nspot.r <- ifposlen(object@maLayout@maNsr)
	MA$printer$nspot.c <- ifposlen(object@maLayout@maNsc)
	MA$printer$notes <- ifposlen(object@maLayout@maNotes)
	MA$genes <- ifposlen(object@maGnames@maInfo)
	MA$genes$Labels <- ifposlen(object@maGnames@maLabels)
	attr(MA$genes,"notes") <- ifposlen(object@maGnames@maNotes)
	MA$genes$Sub <- ifposlen(object@maLayout@maSub)
	MA$genes$Plate <- ifposlen(object@maLayout@maPlate)
	MA$genes$Controls <- ifposlen(object@maLayout@maControls)
	MA$targets <- ifposlen(object@maTargets@maInfo)
	MA$targets$Labels <- ifposlen(object@maTargets@maLabels)
	MA$notes <- ifposlen(object@maNotes)
	MA$maNormCall <- ifposlen(object@maNormCall)
	MA
} 

#  Gordon Smyth, 28 Oct 2004
as.matrix.RGList <- function(x,...) normalizeWithinArrays(x,method="none")$M
as.matrix.MAList <- function(x,...) as.matrix(x$M)
as.matrix.MArrayLM <- function(x,...) x$coefficients
as.matrix.marrayNorm <- function(x,...) x@maM
as.matrix.exprSet <- function(x,...) x@exprs
#  13 July 2006
as.matrix.PLMset <- function(x,...) x@chip.coefs
#  19 Dec 2006, 18 May 2007
as.matrix.ExpressionSet <- as.matrix.LumiBatch <- function(x,...) env=x@assayData[["exprs"]]
#  16 Sep 2007
as.matrix.vsn <- function(x,...) x@hx

if(getRversion() >= "2.4.0") {

as.data.frame.MArrayLM <- function(x, row.names = NULL, optional = FALSE, ...)
#	Convert MAList object to data.frame
#	Gordon Smyth
#	6 April 2005.  Last modified 13 Jan 2006.
{
	x <- unclass(x)
	if(is.null(x$coefficients)) {
		warning("NULL coefficients, returning empty data.frame")
		return(data.frame())
	}
	cn <- names(x)
	nprobes <- NROW(x$coefficients)
	ncoef <- NCOL(x$coefficients)
	include.comp <- cn[unlist(lapply(x,NROW))==nprobes]
	other.comp <- setdiff(names(x),include.comp)
	if(length(other.comp)) for (a in other.comp) x[[a]] <- NULL
#	coef.comp <- c("coefficients","stdev.unscaled","t","p.value","lods")
#	for (a in coef.comp) if(!is.null(x[[a]]) && NCOL(x[[a]])==1) colnames(x[[a]]) <- paste(a,colnames(x[[a]]),sep=".")
	if(ncoef==1) x <- lapply(x,drop)
	as.data.frame(x,row.names=row.names,optional=optional)
}

} else {

as.data.frame.MArrayLM <- function(x, row.names = NULL, optional = FALSE)
#	Convert MAList object to data.frame
#	Gordon Smyth
#	6 April 2005.  Last modified 13 Jan 2006.
{
	x <- unclass(x)
	if(is.null(x$coefficients)) {
		warning("NULL coefficients, returning empty data.frame")
		return(data.frame())
	}
	cn <- names(x)
	nprobes <- NROW(x$coefficients)
	ncoef <- NCOL(x$coefficients)
	include.comp <- cn[unlist(lapply(x,NROW))==nprobes]
	other.comp <- setdiff(names(x),include.comp)
	if(length(other.comp)) for (a in other.comp) x[[a]] <- NULL
#	coef.comp <- c("coefficients","stdev.unscaled","t","p.value","lods")
#	for (a in coef.comp) if(!is.null(x[[a]]) && NCOL(x[[a]])==1) colnames(x[[a]]) <- paste(a,colnames(x[[a]]),sep=".")
	if(ncoef==1) x <- lapply(x,drop)
	as.data.frame(x,row.names=row.names,optional=optional)
}

}
#  MODELMATRIX.R

uniqueTargets <- function(targets) {
	sort(unique(as.vector(as.matrix(targets[,c("Cy3","Cy5")]))))
}

modelMatrix <- function(targets, parameters=NULL, ref=NULL, verbose=TRUE)
#	Design matrix for two-color experiments
#	'targets' is matrix or data.frame with columns Cy3 and Cy5
#	'parameters' specifies desired coefficients corresponding to columns of design matrix
#	'ref' is common reference if such exists
#	Gordon Smyth
#	25 June 2003. Last modified 28 May 2005.
{
	targets <- as.matrix(targets)
	if(missing(targets)) stop("targets is required argument")
	if(!all(c("Cy3","Cy5") %in% colnames(targets))) stop("targets should contain columns: Cy3 and Cy5")
	if(is.null(parameters)==is.null(ref)) stop("exactly one of the arguments parameters and ref should be specified")

	target.names <- sort(unique(as.vector(t(as.matrix(targets[,c("Cy3","Cy5")])))))
	if(verbose) cat("Found unique target names:\n",target.names,"\n")
	if(is.null(parameters)) {
		if(!(ref %in% target.names)) stop(paste("\"",ref,"\" not among the target names found",sep=""))
		other.names <- setdiff(target.names,ref)
		target.names <- c(ref,other.names)
		ntargets <- length(target.names)
		parameters <- rbind(-1,diag(ntargets-1))
		rownames(parameters) <- target.names
		colnames(parameters) <- other.names
	} else {
		parameters <- as.matrix(parameters)
		if(length(target.names) != nrow(parameters)) stop("rows of parameters don't match unique target names")
		if(any(sort(target.names)!=sort(rownames(parameters)))) stop("rownames of parameters don't match unique target names")
		target.names <- rownames(parameters)
		ntargets <- nrow(parameters)
		if(ncol(parameters) != ntargets-1) warning("number of parameters should be one less than number of targets")
	}
	narrays <- nrow(targets)
	J <- matrix(rep(target.names,narrays),ntargets,narrays)
	J <- t((t(J) == targets[,"Cy5"]) - (t(J) == targets[,"Cy3"]))
	rownames(J) <- target.names
	colnames(J) <- rownames(targets)
	zapsmall(t(solve(crossprod(parameters),crossprod(parameters,J))),14)
}

makeContrasts <- function(..., levels)
#	Construct matrix of custom contrasts
#	Gordon Smyth
#	30 June 2003.  Last modified 12 June 2005.
{
	if(is.factor(levels)) levels <- levels(levels)
	if(is.matrix(levels)) levels <- colnames(levels)
	if(is.data.frame(levels)) levels <- names(levels)
	levels <- make.names(levels)
	n <- length(levels)
	if(n < 1) stop("No levels to construct contrasts from")
	indicator <- function(i,n) {
		out <- rep(0,n)
		out[i] <- 1
		out
	}
	for (i in 1:n) assign(levels[i], indicator(i,n))
	e <- substitute(list(...))
	ne <- length(e)
	cm <- matrix(0,n,ne-1)
	rownames(cm) <- levels
	if(ne < 2) return(cm)
	enames <- names(e)[2:ne]
	easchar <- as.character(e)[2:ne]
	if(is.null(enames))
		colnames(cm) <- easchar
	else {
		colnames(cm) <- ifelse(enames=="",easchar,enames)
	}
	for (j in 1:(ne-1)) {
		ej <- e[[j+1]]
		if(is.character(ej)) ej <- parse(text=ej)
		ej <- eval(ej)
#		was original argument a variable?
		if(!is.numeric(ej)) {
			colnames(cm)[j] <- as.character(ej)
			if(is.character(ej)) ej <- parse(text=ej)
			ej <- eval(ej)
		}
		cm[,j] <- ej
	}
	cm
}


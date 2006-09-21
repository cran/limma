#  Quality weights

wtarea <- function(ideal=c(160,170))
#	Quality weights based on spot area from SPOT output
#	Gordon Smyth
#	9 March 2003.  Last revised 11 Mar 2003.

function(spot) {
	e <- range(ideal)
	x <- c(-Inf,0,e,sum(e),Inf)
	y <- c(0,0,1,1,0,0)
	approx(x,y,xout=spot[,"area"],ties="ordered")$y
}

wtflags <- function(weight=0,cutoff=0)
#	Quality weights based on Flags from GenePix output
#	Gordon Smyth
#	9 March 2003.  Last revised 29 July 2006.

function(gpr) {
	flagged <- (gpr[,"Flags"] < cutoff)
	weight*flagged + !flagged
}

wtIgnore.Filter <- function(qta)
#	Quality weights based on Ignore Filter from QuantArray output
#	Gordon Smyth
#	23 May 2003.  Last modified 27 Sep 2003.
{
	qta[,"Ignore Filter"]
}

modifyWeights <- function(weights=rep(1,length(status)), status, values, multipliers)
#	Modify weights for given status values
#	Gordon Smyth
#	29 Dec 2003. Last modified 2 September 2006.
{
	status <- as.character(status)
	weights <- as.matrix(weights)
	values <- as.character(values)
	multipliers <- as.numeric(multipliers)
	if(length(status)!=nrow(weights)) stop("nrows of weights must equal length of status")
	nvalues <- length(values)
	if(length(multipliers)==1) multipliers <- rep(multipliers,nvalues)
	if(nvalues!=length(multipliers)) stop("no. values doesn't match no. multipliers")
	for (i in 1:nvalues) {
		g <- status==values[i]
		weights[g,] <- multipliers[i]*weights[g,]
	}
	weights
}

asMatrixWeights <- function(weights,dim=NULL)
#	Convert probe-weights or array-weights to weight matrix
#	Gordon Smyth
#	22 Jan 2006.
{
	weights <- as.matrix(weights)
	if(is.null(dim)) return(weights)
	if(length(dim)<2) stop("dim must be numeric vector of length 2")
	dim <- round(dim[1:2])
	if(any(dim<1)) stop("zero or negative dimensions not allowed")
	dw <- dim(weights)
#	Full matrix already
	if(all(dw==dim)) return(weights)
	if(min(dw)!=1) stop("weights is of unexpected shape")
#	Row matrix of array weights
	if(dw[2]>1 && dw[2]==dim[2]) return(matrix(weights,dim[1],dim[2],byrow=TRUE))
	lw <- prod(dw)
#	Probe weights
	if(lw==1 || lw==dim[1]) return(matrix(weights,dim[1],dim[2]))
#	Array weights
	if(lw==dim[2]) return(matrix(weights,dim[1],dim[2],byrow=TRUE))
#	All other cases
	stop("weights is of unexpected size")
}



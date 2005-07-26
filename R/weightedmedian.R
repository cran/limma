weighted.median <- function (x, w, na.rm = FALSE)
#	Weighted median
#	Gordon Smyth
#	30 June 2005
{
	if (missing(w)) 
		w <- rep.int(1, length(x))
	else {
		if(length(w) != length(x)) stop("'x' and 'w' must have the same length")
		if(any(is.na(w))) stop("NA weights not allowed")
		if(any(w<0)) stop("Negative weights not allowed")
	}
	if(is.integer(w)) 
		w <- as.numeric(w)
	if(na.rm) {
		w <- w[i <- !is.na(x)]
		x <- x[i]
	}
	if(all(w==0)) {
		warning("All weights are zero")
		return(NA)
	}
	o <- order(x)
	x <- x[o]
	w <- w[o]
	p <- cumsum(w)/sum(w)
	n <- sum(p<0.5)
	if(p[n+1] > 0.5)
		x[n+1]
	else
		(x[n+1]+x[n+2])/2
}

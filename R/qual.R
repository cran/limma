#  QUALITY MEASURES

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

wtflags <- function(w=0.1)
#	Quality weights based on Flags from GenePix output
#	Gordon Smyth
#	9 March 2003.  Last revised 11 June 2003.

function(gpr) {
	flagged <- (gpr[,"Flags"] < 0)
	w*flagged + !flagged
}

wtIgnore.Filter <- function(qta) {
#	Quality weights based on Ignore Filter from QuantArray output
#	Gordon Smyth
#	23 May 2003.  Last modified 27 Sep 2003.

	qta[,"Ignore Filter"]
}

arrayWeightsQuick <- function(y, fit)
#	Compute approximate array quality weights
#	Gordon Smyth
#	25 Oct 2004.  Last revised 28 Oct 2004.
{
	if(!is.null(fit$weights)) warning("spot quality weights found but not taken into account")
	res <- as.matrix(y)- fit$coef %*% t(fit$design)
	h <- hat(fit$design, intercept=FALSE)
	mures2 <- fit$sigma^2 %*% array(1-h,c(1,length(h)))
	1/colMeans(res*res/mures2,na.rm=TRUE)
}



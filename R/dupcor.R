#  DUPCOR.R

duplicateCorrelation <- function(object,design=rep(1,ncol(M)),ndups=2,spacing=1,block=NULL,trim=0.15,weights=NULL)
{
#	Estimate the correlation between duplicates given a series of arrays
#	Gordon Smyth
#	25 Apr 2002. Last revised 5 April 2004.

	if(is(object,"MAList")) {
		M <- object$M
		if(missing(design) && !is.null(object$design)) design <- object$design
		if(missing(ndups) && !is.null(object$printer$ndups)) ndups <- object$printer$ndups
		if(missing(spacing) && !is.null(object$printer$spacing)) spacing <- object$printer$spacing
		if(missing(weights) && !is.null(object$weights)) weights <- object$weights
	} else {
		M <- as.matrix(object)
	}
	narrays <- ncol(M)
	if(is.null(block)) {
		if(ndups<2) {
			warning("No duplicates: correlation between duplicates not estimable")
			return( list(cor=NA,cor.genes=rep(NA,nrow(M))) )
		}
		if(is.character(spacing)) {
			if(spacing=="columns") spacing <- 1
			if(spacing=="rows") spacing <- object$printer$nspot.c
			if(spacing=="topbottom") spacing <- nrow(M)/2
		}
		Array <- rep(1:narrays,rep(ndups,narrays))
	} else {
		ndups <- 1
		nspacing <- 1
		Array <- block
	}
	require( "statmod" ) # need randomizedBlockFit function
	design <- as.matrix(design)
	if(nrow(design) != narrays) stop("Number of rows of design matrix does not match number of arrays")
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		if(any(dim(weights) != dim(M))) weights <- array(weights,dim(M))
		M[weights < 1e-15 ] <- NA
		weights[weights < 1e-15] <- NA
	}
	nbeta <- ncol(design)
	M <- unwrapdups(M,ndups=ndups,spacing=spacing)
	ngenes <- nrow(M)
	if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	design <- design %x% rep(1,ndups)
	rho <- rep(NA,ngenes)
	for (i in 1:ngenes) {
		y <- drop(M[i,])
		o <- is.finite(y)
		A <- Array[o]
		if(sum(o)>(nbeta+2) && any(duplicated(A))) {
			y <- y[o]
			X <- design[o,]
			Z <- model.matrix(~factor(A)-1)
			if(!is.null(weights)) {
				w <- drop(weights[i,])[o]
				s <- randomizedBlockFit(y,X,Z,w,only.varcomp=TRUE,maxit=20)$varcomp
			} else
				s <- randomizedBlockFit(y,X,Z,only.varcomp=TRUE,maxit=20)$varcomp
			rho[i] <- s[2]/sum(s)
		}
	}
	rho <- pmax(-1,rho)
	rhom <- tanh(mean(atanh(rho),trim=trim,na.rm=TRUE))
	list(cor=rhom,consensus.correlation=rhom,all.correlations=rho)
}

dupcor.series <- function(M,design=rep(1,ncol(M)),ndups=2,spacing=1,initial=0.8,trim=0.15,weights=NULL)
{
#	Estimate the correlation between duplicates given a series of arrays
#	Gordon Smyth
#	25 Apr 2002.  Last modified 5 May 2004.
#	This function is deprecated 2 Feb 2004.

	.Deprecated("duplicateCorrelation")
	m <- as.list(match.call())
	m[[1]] <- as.name("duplicateCorrelation")
	m$object <- m$M
	m$M <- m$initial <- NULL
	eval(as.call(m))
}

#  DUPS.R
#  Functions to handle duplicate spots

unwrapdups <- function(M,ndups=2,spacing=1) {
#	Unwrap M matrix for a series of experiments so that all spots for a given gene are in one row
#	Gordon Smyth
#	18 Jan 2002. Last revised 2 Nov 2002.

	if(ndups==1) return(M)
	M <- as.matrix(M)
	nspots <- dim(M)[1]
	nslides <- dim(M)[2]
	ngroups <- nspots / ndups / spacing
	dim(M) <- c(spacing,ndups,ngroups,nslides)
	M <- aperm(M,perm=c(1,3,2,4))
	dim(M) <- c(spacing*ngroups,ndups*nslides)
	M
}

uniquegenelist <- function(genelist,ndups=2,spacing=1) {
#	Eliminate entries in genelist for duplicate spots
#	Gordon Smyth
#	2 Nov 2002.  Last revised 10 Jan 2005

	if(ndups <= 1) return(genelist)
	i <- drop(unwrapdups(1:NROW(genelist),ndups=ndups,spacing=spacing)[,1])
	if(is.null(dim(genelist)))
		return(genelist[i])
	else
		return(genelist[i,,drop=FALSE])
}

duplicateCorrelation <- function(object,design=rep(1,ncol(as.matrix(object))),ndups=2,spacing=1,block=NULL,trim=0.15,weights=NULL)
#	Estimate the correlation between duplicates given a series of arrays
#	Gordon Smyth
#	25 Apr 2002. Last revised 14 Nov 2005.
{
	M <- as.matrix(object)
	if(is(object,"MAList")) {
		if(missing(design) && !is.null(object$design)) design <- object$design
		if(missing(ndups) && !is.null(object$printer$ndups)) ndups <- object$printer$ndups
		if(missing(spacing) && !is.null(object$printer$spacing)) spacing <- object$printer$spacing
		if(missing(weights) && !is.null(object$weights)) weights <- object$weights
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
	require( "statmod" ) # need mixedModel2Fit function
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
		A <- factor(Array[o])
		nobs <- sum(o)
		nblocks <- length(levels(A))
		if(nobs>(nbeta+2) && nblocks>1 && nblocks<nobs-1) {
			y <- y[o]
			X <- design[o,]
			Z <- model.matrix(~0+A)
			if(!is.null(weights)) {
				w <- drop(weights[i,])[o]
				s <- mixedModel2Fit(y,X,Z,w,only.varcomp=TRUE,maxit=20)$varcomp
			} else
				s <- mixedModel2Fit(y,X,Z,only.varcomp=TRUE,maxit=20)$varcomp
			rho[i] <- s[2]/sum(s)
		}
	}
	arho <- atanh(pmax(-1,rho))
	mrho <- tanh(mean(arho,trim=trim,na.rm=TRUE))
	list(consensus.correlation=mrho,cor=mrho,atanh.correlations=arho)
}

avedups <- function(x,ndups,spacing,weights) UseMethod("avedups")

avedups.default <- function(x,ndups=2,spacing=1,weights=NULL)
#	Average over duplicate spots, for matrices or vectors
#	Gordon Smyth
#	6 Apr 2006.
{
	if(ndups==1) return(x)
	if(is.null(x)) return(NULL)
	x <- as.matrix(x)
	nspots <- dim(x)[1]
	nslides <- dim(x)[2]
	rn <- rownames(x)
	cn <- colnames(x)
	ngroups <- nspots / ndups / spacing
	dim(x) <- c(spacing,ndups,ngroups*nslides)
	x <- aperm(x,perm=c(2,1,3))
	if(mode(x)=="character")
		x <- x[1,,]
	else {
		if(is.null(weights))
			x <- colMeans(x,na.rm=TRUE)
		else {
			weights <- as.matrix(weights)
			dim(weights) <- c(spacing,ndups,ngroups*nslides)
			weights <- aperm(weights,perm=c(2,1,3))
			weights[is.na(weights) | is.na(x)] <- 0
			weights[weights<0] <- 0
			x <- colSums(weights*x,na.rm=TRUE)/colSums(weights)
		}
	}
	dim(x) <- c(spacing*ngroups,nslides)
	colnames(x) <- cn
	rownames(x) <- avedups(rn,ndups=ndups,spacing=spacing)
	x
}

avedups.MAList <- function(x,ndups=x$printer$ndups,spacing=x$printer$spacing,weights=x$weights)
#	Average over duplicate spots for MAList objects
#	Gordon Smyth
#	6 Apr 2006.
{
	if(is.null(ndups) || is.null(spacing)) stop("Must specify ndups and spacing")
	y <- x
	y$M <- avedups(x$M,ndups=ndups,spacing=spacing,weights=weights)
	y$A <- avedups(x$A,ndups=ndups,spacing=spacing,weights=weights)
	other <- names(x$other)
	for (a in other) object$other[[a]] <- avedups(object$other[[a]],ndups=ndups,spacing=spacing,weights=weights)
	y$weights <- avedups(x$weights,ndups=ndups,spacing=spacing)
	y$genes <- uniquegenelist(x$genes,ndups=ndups,spacing=spacing)
    y$printer <- NULL
	y
}


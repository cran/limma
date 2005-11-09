#  LINEAR MODELS

lmFit <- function(object,design=NULL,ndups=1,spacing=1,block=NULL,correlation,weights=NULL,method="ls",...) {
#	Fit linear model
#	Gordon Smyth
#	30 June 2003.  Last modified 7 Nov 2005.

	M <- NULL
#	Method intended for MAList objects but allow unclassed lists as well
	if(is(object,"MAList") || is(object,"list")) {
		M <- object$M
		if(missing(design) && !is.null(object$design)) design <- object$design
		if(missing(ndups) && !is.null(object$printer$ndups)) ndups <- object$printer$ndups
		if(missing(spacing) && !is.null(object$printer$spacing)) spacing <- object$printer$spacing
		if(missing(correlation) && !is.null(object$correlation)) correlation <- object$correlation
		if(missing(weights) && !is.null(object$weights)) weights <- object$weights
	} else {
	if(is(object,"marrayNorm")) {
#		don't use accessor function so don't have to require marrayClasses
		M <- object@maM
		if(missing(weights) && length(object@maW)) weights <- object@maW
	} else {
	if(is(object,"PLMset")) {
#		don't use accessor function so don't have to require affyPLM
		M <- object@chip.coefs
		if(length(M)==0) stop("chip.coefs has length zero")
		if(missing(weights) && length(object@se.chip.coefs)) weights <- 1/pmax(object@se.chip.coefs,1e-5)^2
	} else {
	if(is(object,"exprSet")) {
#		don't use accessor function so don't have to require Biobase
		M <- object@exprs
#		don't use weights until this is more thoroughly tested
#		if(missing(weights) && length(object@se.exprs)) weights <- 1/pmax(object@se.exprs,1e-5)^2
	}}}}
#	Default method
	if(is.null(M)) M <- as.matrix(object)

	if(is.null(design)) design <- matrix(1,ncol(M),1)
	design <- as.matrix(design)
	ne <- nonEstimable(design)
	if(!is.null(ne)) cat("Coefficients not estimable:",paste(ne,collapse=" "),"\n")
	method <- match.arg(method,c("ls","robust"))
	if(method=="robust")
		fit <- mrlm(M,design=design,ndups=ndups,spacing=spacing,weights=weights,...)
	else
		if(ndups < 2 && is.null(block))
			fit <- lm.series(M,design=design,ndups=ndups,spacing=spacing,weights=weights)
		else {
			if(missing(correlation)) stop("the correlation must be set, see duplicateCorrelation")
			fit <- gls.series(M,design=design,ndups=ndups,spacing=spacing,block=block,correlation=correlation,weights=weights,...)
		}
	fit$method <- method
	fit$design <- design
	if(is(object,"MAList")) {
		if(!is.null(object$genes)) fit$genes <- uniquegenelist(object$genes,ndups=ndups,spacing=spacing) 
		if(!is.null(object$A)) fit$Amean <- rowMeans(unwrapdups(as.matrix(object$A),ndups=ndups,spacing=spacing),na.rm=TRUE)
	}
	if(is(object,"marrayNorm")) {
		if(length(object@maGnames@maInfo)) {
			fit$genes <- object@maGnames@maInfo
#			if(length(object@maLayout@maControls)>1) fit$genes$Status <- normdata@maLayout@maControls
			fit$genes <- uniquegenelist(fit$genes,ndups=ndups,spacing=spacing)
			attr(fit$genes, "Notes") <- object@maGnames@maNotes
		}
		if(length(object@maA)) fit$Amean <- rowMeans(unwrapdups(object@maA,ndups=ndups,spacing=spacing),na.rm=TRUE)
	}
	if(is(object,"exprSet") || is(object,"matrix")) {
		ProbeID <- rownames(M)
		if(!is.null(ProbeID)) fit$genes <- uniquegenelist(data.frame(ID=I(ProbeID)),ndups=ndups,spacing=spacing)
		fit$Amean <- rowMeans(M,na.rm=TRUE)
	}
	if(is.null(fit$genes)) fit$genes <- rownames(M)
	new("MArrayLM",fit)
}

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

lm.series <- function(M,design=NULL,ndups=1,spacing=1,weights=NULL)
{
#	Fit linear model for each gene to a series of arrays
#	Gordon Smyth
#	18 Apr 2002. Revised 13 January 2005.

	M <- as.matrix(M)
	narrays <- ncol(M)
	if(is.null(design)) design <- matrix(1,narrays,1)
	design <- as.matrix(design)
	nbeta <- ncol(design)
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		if(any(dim(weights) != dim(M))) weights <- array(weights,dim(M))
		weights[weights <= 0] <- NA
		M[!is.finite(weights)] <- NA
	}
	if(ndups>1) {
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		design <- design %x% rep(1,ndups)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	}
	ngenes <- nrow(M)
	stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta,dimnames=list(rownames(M),colnames(design)))
	sigma <- rep(NA,ngenes)
	df.residual <- rep(0,ngenes)
	NoWts <- !any(is.na(M)) && is.null(weights)
	if(NoWts) {
		fit <- lm.fit(design, t(M))
		if(fit$df.residual>0)
			fit$sigma <- sqrt(colMeans(fit$effects[(fit$rank + 1):narrays,,drop=FALSE]^2))
		else
			fit$sigma <- rep(NA,ngenes)
		fit$fitted.values <- fit$residuals <- fit$effects <- NULL
		fit$coefficients <- t(fit$coefficients)
		fit$cov.coefficients <- chol2inv(fit$qr$qr,size=fit$qr$rank)
		est <- fit$qr$pivot[1:fit$qr$rank]
		stdev.unscaled[,est] <- matrix(sqrt(diag(fit$cov.coefficients)),ngenes,fit$qr$rank,byrow = TRUE)
		fit$stdev.unscaled <- stdev.unscaled
		fit$df.residual <- rep.int(fit$df.residual,ngenes)
		dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
		fit$pivot <- fit$qr$pivot
		return(fit)
	}
	for (i in 1:ngenes) {
		y <- as.vector(M[i,])
		obs <- is.finite(y)
		if(sum(obs) > 0) {
			X <- design[obs,,drop=FALSE]
			y <- y[obs]
			if(is.null(weights))
				out <- lm.fit(X,y)
			else {
				w <- as.vector(weights[i,obs])
				out <- lm.wfit(X,y,w)
			}
			est <- !is.na(out$coef)
			beta[i,] <- out$coef
			stdev.unscaled[i,est] <- sqrt(diag(chol2inv(out$qr$qr,size=out$rank)))
			df.residual[i] <- out$df.residual
			if(df.residual[i] > 0)
				if(is.null(weights))
					sigma[i] <- sqrt(sum(out$residuals^2)/out$df.residual)
				else
					sigma[i] <- sqrt(sum(w*out$residuals^2)/out$df.residual)
		}
	}
	QR <- qr(design)
	cov.coef <- chol2inv(QR$qr,size=QR$rank)
	list(coefficients=drop(beta),stdev.unscaled=drop(stdev.unscaled),sigma=sigma,df.residual=df.residual,cov.coefficients=cov.coef,pivot=QR$pivot)
}

rlm.series <- function(x,...)
#	Gordon Smyth
#	Deprecated 1 Sep 2004
{
	.Deprecated("mrlm")
	m <- match.call()
	m[[1]] <- as.name("mrlm")
	names(m)[2] <- "M"
	eval(m)
}

mrlm <- function(M,design=NULL,ndups=1,spacing=1,weights=NULL,...)
{
#	Robustly fit linear model for each gene to a series of arrays
#	Gordon Smyth
#	20 Mar 2002.  Last revised 19 Sept 2003.

	require(MASS) # need rlm.default
	M <- as.matrix(M)
	narrays <- ncol(M)
	if(is.null(design)) design <- matrix(1,narrays,1)
	design <- as.matrix(design)
	nbeta <- ncol(design)
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		if(any(dim(weights) != dim(M))) weights <- array(weights,dim(M))
		weights[weights <= 0] <- NA
		M[!is.finite(weights)] <- NA
	}
	if(ndups>1) {
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		design <- design %x% rep(1,ndups)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	}
	ngenes <- nrow(M)
	stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta)
	sigma <- rep(NA,ngenes)
	df.residual <- rep(0,ngenes)
	for (i in 1:ngenes) {
		y <- as.vector(M[i,])
		obs <- is.finite(y)
		X <- design[obs,,drop=FALSE]
		y <- y[obs]
		if(is.null(weights))
			w <- rep(1,length(y))
		else
			w <- as.vector(weights[i,obs])
		if(length(y) > nbeta) {
			out <- rlm(x=X,y=y,weights=w,...)
			beta[i,] <- coef(out)
			stdev.unscaled[i,] <- sqrt(diag(chol2inv(out$qr$qr)))
			df.residual[i] <- length(y) - out$rank
			if(df.residual[i] > 0) sigma[i] <- out$s
		}
	}
	QR <- qr(design)
	cov.coef <- chol2inv(QR$qr,size=QR$rank)
	list(coefficients=drop(beta),stdev.unscaled=drop(stdev.unscaled),sigma=sigma,df.residual=df.residual,cov.coefficients=cov.coef,pivot=QR$pivot)
}

gls.series <- function(M,design=NULL,ndups=2,spacing=1,block=NULL,correlation=NULL,weights=NULL,...)
{
#	Fit linear model for each gene to a series of microarrays.
#	Fit is by generalized least squares allowing for correlation between duplicate spots.
#	Gordon Smyth
#	11 May 2002.  Last revised 21 April 2005.

	M <- as.matrix(M)
	narrays <- ncol(M)
	if(is.null(design)) design <- matrix(1,narrays,1)
	design <- as.matrix(design)
	if(nrow(design) != narrays) stop("Number of rows of design matrix does not match number of arrays")
	if(is.null(correlation)) correlation <- duplicateCorrelation(M,design=design,ndups=ndups,spacing=spacing,block=block,weights=weights,...)$consensus.correlation
	if(is.null(block)) {
		if(ndups<2) {
			warning("No duplicates: correlation between duplicates set to zero")
			ndups <- 1
			correlation <- 0
		}
		if(is.null(spacing)) spacing <- 1
		cormatrix <- diag(rep(correlation,len=narrays)) %x% array(1,c(ndups,ndups))
	} else {
		if(ndups>1) {
			stop("Cannot specify ndups>2 and non-null block argument")
		} else {
			ndups <- spacing <- 1
		}
		block <- as.vector(block)
		if(length(block)!=narrays) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,narrays,nblocks)==matrix(ub,narrays,nblocks,byrow=TRUE)
		cormatrix <- Z%*%(correlation*t(Z))
	}
	diag(cormatrix) <- 1
	if(!is.null(weights)) {
		weights <- as.matrix(weights)
		if(any(dim(weights) != dim(M))) weights <- array(weights,dim(M))
		M[weights < 1e-15 ] <- NA
		weights[weights < 1e-15] <- NA
	}
	nbeta <- ncol(design)
	coef.names <- colnames(design)
	M <- unwrapdups(M,ndups=ndups,spacing=spacing)
	ngenes <- nrow(M)
	if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	design <- design %x% rep(1,ndups)
	stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta,dimnames=list(NULL,coef.names))
	sigma <- rep(NA,ngenes)
	df.residual <- rep(0,ngenes)
	for (i in 1:ngenes) {
		y <- drop(M[i,])
		o <- is.finite(y)
		y <- y[o]
		n <- length(y)
		if(n > 0) {
			X <- design[o,,drop=FALSE]
			V <- cormatrix[o,o]
			if(!is.null(weights)) {
				wrs <- 1/sqrt(drop(weights[i,o]))
				V <- wrs * t(wrs * t(V))
			}
			cholV <- chol(V)
			y <- backsolve(cholV,y,transpose=TRUE)
			if(all(X==0)) {
				df.residual[i] <- n
				sigma[i] <- sqrt( array(1/n,c(1,n)) %*% y^2 )
			} else {
				X <- backsolve(cholV,X,transpose=TRUE)
				out <- lm.fit(X,y)
				est <- !is.na(out$coefficients)
				beta[i,] <- out$coefficients
				stdev.unscaled[i,est] <- sqrt(diag(chol2inv(out$qr$qr,size=out$rank)))
				df.residual[i] <- out$df.residual
				if(df.residual[i] > 0)
					sigma[i] <- sqrt( array(1/out$df.residual,c(1,n)) %*% out$residuals^2 )
			}
		}
	}
	cholV <- chol(cormatrix)
	QR <- qr(backsolve(cholV,design,transpose=TRUE))
	cov.coef <- chol2inv(QR$qr,size=QR$rank)
	list(coefficients=drop(beta),stdev.unscaled=drop(stdev.unscaled),sigma=sigma,df.residual=df.residual,ndups=ndups,spacing=spacing,block=block,correlation=correlation,cov.coefficients=cov.coef,pivot=QR$pivot)
}

contrasts.fit <- function(fit,contrasts) {
#	Convert coefficients and std deviations in fit object to reflect contrasts of interest
#	Gordon Smyth
#	13 Oct 2002.  Last modified 25 Sep 2004.

	ncoef <- NCOL(fit$coefficients)
	if(NROW(contrasts)!=ncoef) stop("Number of rows of contrast matrix must match number of coefficients")
	fit$contrasts <- contrasts
	cormatrix <- cov2cor(fit$cov.coefficients)
	if(is.null(cormatrix)) {
		warning("no coef correlation matrix found in fit - assuming orthogonal")
		cormatrix <- diag(ncoef)
	}

#	If design matrix was singular, reduce to estimable coefficients
	r <- nrow(cormatrix)
	if(r < ncoef) {
		if(is.null(fit$pivot)) stop("cor.coef not full rank but pivot column not found in fit")
		est <- fit$pivot[1:r]
		if(any(contrasts[-est,])) stop("trying to take contrast of non-estimable coefficient")
		contrasts <- contrasts[est,,drop=FALSE]
		fit$coefficients <- fit$coefficients[,est,drop=FALSE]
		fit$stdev.unscaled <- fit$stdev.unscaled[,est,drop=FALSE]
		ncoef <- r
	}
	fit$coefficients <- fit$coefficients %*% contrasts

#	Test whether design was orthogonal
	if(length(cormatrix) < 2) {
		orthog <- TRUE
	} else {
		orthog <- all(abs(cormatrix[lower.tri(cormatrix)]) < 1e-14)
	}

#	Contrast correlation matrix
	R <- La.chol(fit$cov.coefficients)
	fit$cov.coefficients <- crossprod(R %*% contrasts)
	fit$pivot <- NULL

	if(orthog)
		fit$stdev.unscaled <- sqrt(fit$stdev.unscaled^2 %*% contrasts^2)
	else {
		R <- La.chol(cormatrix)
		ngenes <- NROW(fit$stdev.unscaled)
		ncont <- NCOL(contrasts)
		U <- matrix(1,ngenes,ncont,dimnames=list(rownames(fit$stdev.unscaled),colnames(contrasts)))
		o <- array(1,c(1,ncoef))
		for (i in 1:ngenes) {
			RUC <- R %*% vecmat(fit$stdev.unscaled[i,],contrasts)
			U[i,] <- sqrt(o %*% RUC^2)
		}
		fit$stdev.unscaled <- U
	}
	fit
}

#contrasts.fit0 <- function(fit,contrasts,design=NULL) {
##	Convert coefficients and std deviations in fit object to reflect contrasts of interest
##	Gordon Smyth
##	13 Oct 2002.  Last modified 20 May 2004.
#
#	ncoef <- NCOL(fit$coefficients)
#	if(nrow(contrasts)!=ncoef) stop("Number of rows of contrast matrix must match number of coefficients")
#	fit$coefficients <- fit$coefficients %*% contrasts
#	if(is.null(design)) design <- fit$design
#	if(!is.null(design) && ncoef > 1) {
#		A <- crossprod( abs(design) > 1e-14 )
#		orthog <- all(A[lower.tri(A)]==0) 
#	}
#	if(is.null(design) || ncoef==1 || orthog)
#		fit$stdev.unscaled <- sqrt(fit$stdev.unscaled^2 %*% contrasts^2)
#	else {
#		A <- La.chol2inv(La.chol(crossprod(design)))
#		s <- sqrt(diag(A))
#		R <- La.chol(t(A/s)/s)
#		ngenes <- NROW(fit$stdev.unscaled)
#		ncont <- NCOL(contrasts)
#		U <- matrix(1,ngenes,ncont,dimnames=list(rownames(fit$stdev.unscaled),colnames(contrasts)))
#		o <- array(1,c(1,ncoef))
#		for (i in 1:ngenes) {
#			RUC <- R %*% vecmat(fit$stdev.unscaled[i,],contrasts)
#			U[i,] <- sqrt(o %*% RUC^2)
#		}
#		fit$stdev.unscaled <- U
#	}
#	fit$contrasts <- contrasts
#	fit
#}

#contrasts.fit <- function(fit,contrasts) {
##	Extract contrast information from oneway linear model fit
##	Gordon Smyth
##	13 Oct 2002.  Last modified 1 July 2003.
#
#	fit$coefficients <- fit$coefficients %*% contrasts
#	fit$stdev.unscaled <- sqrt(fit$stdev.unscaled^2 %*% contrasts^2)
#	fit$contrasts <- contrasts
#	fit
#}

is.fullrank <- function(x)
#	Check whether a numeric matrix has full column rank
#	Gordon Smyth
#	18 August 2003.  Last modified 9 March 2004.
{
	x <- as.matrix(x)
	e <- eigen(crossprod(x),symmetric=TRUE,only.values=TRUE)$values
	e[1] > 0 && abs(e[length(e)]/e[1]) > 1e-13
}

nonEstimable <- function(x)
#	Check whether a numeric matrix has full column rank
#	If not, return names of redundant columns
#	Gordon Smyth
#	10 August 2004
{
	x <- as.matrix(x)
	p <- ncol(x)
	QR <- qr(x)
	if(QR$rank < p) {
		n <- colnames(x)
		if(is.null(n)) n <- as.character(1:p)
		notest <- n[QR$pivot[(QR$rank+1):p]]
		blank <- notest==""
		if(any(blank)) notest[blank] <- as.character(((QR$rank+1):p)[blank])
		return(notest)
	} else {
		return(NULL)
	}
}

#  LINEAR MODELS

lmFit <- function(object,design=NULL,ndups=1,spacing=1,block=NULL,correlation,weights=NULL,method="ls",...)
#	Fit linear model
#	Gordon Smyth
#	30 June 2003.  Last modified 10 March 2007.
{
	y <- NULL
	Amean <- NULL
	ProbeAnn <- NULL
	if(is(object,"list")) {
#		Method intended for MAList objects but allow unclassed lists as well
		y <- object$M
		if(missing(design) && !is.null(object$design)) design <- object$design
		if(missing(ndups) && !is.null(object$printer$ndups)) ndups <- object$printer$ndups
		if(missing(spacing) && !is.null(object$printer$spacing)) spacing <- object$printer$spacing
		if(missing(correlation) && !is.null(object$correlation)) correlation <- object$correlation
		if(missing(weights) && !is.null(object$weights)) weights <- object$weights
		ProbeAnn <- object$genes
		if(!is.null(object$A)) Amean <- rowMeans(unwrapdups(as.matrix(object$A),ndups=ndups,spacing=spacing),na.rm=TRUE)
	} else {
	if(is(object,"marrayNorm")) {
#		don't use accessor function so don't have to require marrayClasses
		y <- object@maM
		if(missing(weights) && length(object@maW)) weights <- object@maW
		if(length(object@maGnames@maInfo)) {
			ProbeAnn <- object@maGnames@maInfo
#			if(length(object@maLayout@maControls)>1) ProbeAnn$Status <- normdata@maLayout@maControls
			attr(ProbeAnn, "Notes") <- object@maGnames@maNotes
		}
		if(length(object@maA)) Amean <- rowMeans(unwrapdups(object@maA,ndups=ndups,spacing=spacing),na.rm=TRUE)
	} else {
	if(is(object,"PLMset")) {
#		don't use accessor function so don't have to require affyPLM
		y <- object@chip.coefs
		if(length(y)==0) stop("chip.coefs has length zero")
		if(!is.null(rownames(y))) ProbeAnn <- data.frame(ID=I(rownames(y)))
		if(missing(weights) && length(object@se.chip.coefs)) weights <- 1/pmax(object@se.chip.coefs,1e-5)^2
		Amean <- rowMeans(unwrapdups(y,ndups=ndups,spacing=spacing),na.rm=TRUE)
	} else {
	if(is(object,"exprSet")) {
#		don't use accessor function so don't have to require Biobase
		y <- object@exprs
		if(!is.null(rownames(y))) ProbeAnn <- data.frame(ID=I(rownames(y)))
#		don't use weights until this is more thoroughly tested
#		if(missing(weights) && length(object@se.exprs)) weights <- 1/pmax(object@se.exprs,1e-5)^2
		Amean <- rowMeans(unwrapdups(y,ndups=ndups,spacing=spacing),na.rm=TRUE)
	} else {
	if(is(object,"ExpressionSet")) {
		y <- get("exprs",env=object@assayData)
		if(length(object@featureData@data)) ProbeAnn <- object@featureData@data
		if(!is.null(rownames(y))) {
			if(is.null(ProbeAnn))
				ProbeAnn <- data.frame(ID=I(rownames(y)))
			else
				ProbeAnn$ID <- rownames(y)
		}
		Amean <- rowMeans(unwrapdups(y,ndups=ndups,spacing=spacing),na.rm=TRUE)
	}}}}}
#	Default method
	if(is.null(y)) {
		y <- as.matrix(object)
		if(!is.null(rownames(y))) ProbeAnn <- data.frame(ID=I(rownames(y)))
	}

	if(is.null(design)) design <- matrix(1,ncol(y),1)
	design <- as.matrix(design)
	if(mode(design) != "numeric") stop("design must be a numeric matrix")
	ne <- nonEstimable(design)
	if(!is.null(ne)) cat("Coefficients not estimable:",paste(ne,collapse=" "),"\n")

	method <- match.arg(method,c("ls","robust"))
	if(method=="robust")
		fit <- mrlm(y,design=design,ndups=ndups,spacing=spacing,weights=weights,...)
	else
		if(ndups < 2 && is.null(block))
			fit <- lm.series(y,design=design,ndups=ndups,spacing=spacing,weights=weights)
		else {
			if(missing(correlation)) stop("the correlation must be set, see duplicateCorrelation")
			fit <- gls.series(y,design=design,ndups=ndups,spacing=spacing,block=block,correlation=correlation,weights=weights,...)
		}

	if(!is.null(ProbeAnn)) fit$genes <- uniquegenelist(ProbeAnn,ndups=ndups,spacing=spacing)
	fit$Amean <- Amean
	fit$method <- method
	fit$design <- design
	new("MArrayLM",fit)
}

lm.series <- function(M,design=NULL,ndups=1,spacing=1,weights=NULL)
{
#	Fit linear model for each gene to a series of arrays
#	Gordon Smyth
#	18 Apr 2002. Revised 22 Feb 2007.

	M <- as.matrix(M)
	narrays <- ncol(M)
	if(is.null(design)) design <- matrix(1,narrays,1)
	design <- as.matrix(design)
	nbeta <- ncol(design)
	coef.names <- colnames(design)
	if(!is.null(weights)) {
		weights <- asMatrixWeights(weights,dim(M))
		weights[weights <= 0] <- NA
		M[!is.finite(weights)] <- NA
	}
	if(ndups>1) {
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		design <- design %x% rep(1,ndups)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	}
	ngenes <- nrow(M)
	stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta,dimnames=list(rownames(M),coef.names))
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
		dimnames(fit$cov.coefficients) <- list(coef.names[est],coef.names[est])
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
	est <- QR$pivot[1:QR$rank]
	dimnames(cov.coef) <- list(coef.names[est],coef.names[est])
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
#	Robustly fit linear model for each gene to a series of arrays
#	Gordon Smyth
#	20 Mar 2002.  Last revised 22 Feb 2007.
{
	require(MASS) # need rlm.default
	M <- as.matrix(M)
	narrays <- ncol(M)
	if(is.null(design)) design <- matrix(1,narrays,1)
	design <- as.matrix(design)
	coef.names <- colnames(design)
	nbeta <- ncol(design)
	if(!is.null(weights)) {
		weights <- asMatrixWeights(weights,dim(M))
		weights[weights <= 0] <- NA
		M[!is.finite(weights)] <- NA
	}
	if(ndups>1) {
		M <- unwrapdups(M,ndups=ndups,spacing=spacing)
		design <- design %x% rep(1,ndups)
		if(!is.null(weights)) weights <- unwrapdups(weights,ndups=ndups,spacing=spacing)
	}
	ngenes <- nrow(M)
	stdev.unscaled <- beta <- matrix(NA,ngenes,nbeta,dimnames=list(rownames(M),coef.names))
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
	est <- QR$pivot[1:QR$rank]
	dimnames(cov.coef) <- list(coef.names[est],coef.names[est])
	list(coefficients=drop(beta),stdev.unscaled=drop(stdev.unscaled),sigma=sigma,df.residual=df.residual,cov.coefficients=cov.coef,pivot=QR$pivot)
}

gls.series <- function(M,design=NULL,ndups=2,spacing=1,block=NULL,correlation=NULL,weights=NULL,...)
#	Fit linear model for each gene to a series of microarrays.
#	Fit is by generalized least squares allowing for correlation between duplicate spots.
#	Gordon Smyth
#	11 May 2002.  Last revised 22 Feb 2007.
{
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
		weights <- asMatrixWeights(weights,dim(M))
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
	est <- QR$pivot[1:QR$rank]
	dimnames(cov.coef) <- list(coef.names[est],coef.names[est])
	list(coefficients=drop(beta),stdev.unscaled=drop(stdev.unscaled),sigma=sigma,df.residual=df.residual,ndups=ndups,spacing=spacing,block=block,correlation=correlation,cov.coefficients=cov.coef,pivot=QR$pivot)
}

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

fitted.MArrayLM <- function(object,design=object$design,...)
#	Fitted values from MArray linear model fit
#	Gordon Smyth
#	29 November 2005
{
	object$coefficients %*% t(object$design)
}

residuals.MArrayLM <- function(object,y,...)
#	Residuals from MArray linear model fit
#	Gordon Smyth
#	29 November 2005
{
	as.matrix(y) - fitted(object)
}


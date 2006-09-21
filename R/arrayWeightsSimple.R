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

arrayWeightsSimple <- function(object,design=NULL,maxiter=100,tol=1e-6,maxratio=100,trace=FALSE)
#	Array weights by REML
#	Assumes no spot weights
#	Any probes with missing values are removed
#	Gordon Smyth, 13 Dec 2005
#	Last revised 22 Jan 2006.
{
	M <- as.matrix(object)
	allfin <- apply(is.finite(M),1,all)
	if(!all(allfin)) {
		nrowna <- sum(!allfin)
		M <- M[allfin,,drop=FALSE]
		warning(paste(as.integer(nrowna),"rows with missing values removed"))
	}
	ngenes <- nrow(M)
	narrays <- ncol(M)
	if(narrays < 3) stop("too few arrays")
	if(ngenes < narrays) stop("too few probes")
	if(is.null(design)) design <- matrix(1,narrays,1)
	p <- ncol(design)

	Z1 <- contr.sum(narrays)
	Z <- cbind(1,Z1)

#	Starting values
	gam <- rep(0,narrays-1)
	w <- drop(exp(Z1 %*% (-gam)))

	iter <- 0
	p2 <- p*(p+1)/2
	Q2 <- array(0,c(narrays,p2))
	repeat {
		iter <- iter+1
		fitm <- lm.wfit(design, t(M), w)
		Q <- qr.qy(fitm$qr,diag(1,nrow=narrays,ncol=p))
		j0 <- 0
		for (k in 0:(p-1)) {
			Q2[,(j0+1):(j0+p-k)] <- Q[,1:(p-k)]*Q[,(k + 1):p]
			j0 <- j0 + p - k
		}
		if(p > 1) Q2[,(p+1):p2] <- sqrt(2)*Q2[,(p+1):p2]
		h <- rowSums(Q2[,1:p,drop=FALSE])
		info <- crossprod(Z,(1-2*h)*Z) + crossprod(crossprod(Q2,Z))
		info1 <- info[-1,-1,drop=FALSE] - (info[-1,1,drop=FALSE]/info[1,1]) %*% info[1,-1,drop=FALSE]
		s2 <- colSums(fitm$effects[(fitm$rank+1):narrays,]^2) / fitm$df.residual
		fitvald <- matrix(1/w,narrays,1)%*%s2
		dl1 <- crossprod(Z1, rowMeans(fitm$residuals^2/fitvald - (1-h)) )
#		print(cbind(info1,dl1))
#		cat(iter,drop(crossprod(ngenes*dl1)),"\n")
		gamstep <- solve(info1,dl1)
		gam <- gam + gamstep
		w <- drop(exp(Z1 %*% (-gam)))
		convcrit <- ngenes/narrays*crossprod(dl1,gamstep)
		if(trace) cat(iter,convcrit,w,"\n")
		if(convcrit < tol) break
		if(max(w)/min(w) > maxratio) break
		if(iter==maxiter) {
			warning("iteration limit reached")
			break
		}
	}
	w
}

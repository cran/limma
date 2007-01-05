#  CONTRASTS

contrasts.fit <- function(fit,contrasts=NULL,coefficients=NULL) {
#	Convert coefficients and std deviations in fit object to reflect contrasts of interest
#	Gordon Smyth
#	13 Oct 2002.  Last modified 27 August 2006.

	ncoef <- NCOL(fit$coefficients)
	if(is.null(contrasts) == is.null(coefficients)) stop("Must specify only one of contrasts or coefficients")
	if(!is.null(coefficients)) {
		ncont <- length(coefficients)
		contrasts <- diag(ncoef)
		rownames(contrasts) <- colnames(contrasts) <- colnames(fit$coefficients)
		contrasts <- contrasts[,coefficients,drop=FALSE]
	}
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
	R <- chol(fit$cov.coefficients)
	fit$cov.coefficients <- crossprod(R %*% contrasts)
	fit$pivot <- NULL

	if(orthog)
		fit$stdev.unscaled <- sqrt(fit$stdev.unscaled^2 %*% contrasts^2)
	else {
		R <- chol(cormatrix)
		ngenes <- NROW(fit$stdev.unscaled)
		ncont <- NCOL(contrasts)
		U <- matrix(1,ngenes,ncont,dimnames=list(rownames(fit$stdev.unscaled),colnames(contrasts)))
		o <- array(1,c(1,ncoef))
		for (i in 1:ngenes) {
			RUC <- R %*% .vecmat(fit$stdev.unscaled[i,],contrasts)
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
#		A <- chol2inv(chol(crossprod(design)))
#		s <- sqrt(diag(A))
#		R <- chol(t(A/s)/s)
#		ngenes <- NROW(fit$stdev.unscaled)
#		ncont <- NCOL(contrasts)
#		U <- matrix(1,ngenes,ncont,dimnames=list(rownames(fit$stdev.unscaled),colnames(contrasts)))
#		o <- array(1,c(1,ncoef))
#		for (i in 1:ngenes) {
#			RUC <- R %*% .vecmat(fit$stdev.unscaled[i,],contrasts)
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


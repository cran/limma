#  DECIDETESTS.R

setClass("TestResults",representation("matrix"))

summary.TestResults <- function(object,...)
#	Gordon Smyth
#	26 Feb 2004.  Last modified 31 Jan 2005.
{
#	apply(object,2,table)
	tab <- array(0,c(3,ncol(object)),dimnames=list(c("-1","0","1"),colnames(object)))
	tab[1,] <- colSums(object== -1,na.rm=TRUE)
	tab[2,] <- colSums(object== 0,na.rm=TRUE)
	tab[3,] <- colSums(object== 1,na.rm=TRUE)
	class(tab) <- "table"
	tab
}

setMethod("show","TestResults",function(object) {
	cat("TestResults matrix\n")
	printHead(object@.Data)
})

decideTests <- function(object,method="separate",adjust.method="BH",p.value=0.05,lfc=0)
#	Accept or reject hypothesis tests across genes and contrasts
#	Gordon Smyth
#	17 Aug 2004. Last modified 24 June 2007.
{
	if(!is(object,"MArrayLM")) stop("Need MArrayLM object")
	if(is.null(object$t)) object <- eBayes(object)
	method <- match.arg(method,c("separate","global","hierarchical","nestedF"))
	adjust.method <- match.arg(adjust.method,c("none","bonferroni","holm","BH","fdr","BY"))
	if(adjust.method=="fdr") adjust.method <- "BH"
	switch(method,separate={
		p <- as.matrix(object$p.value)
		tstat <- as.matrix(object$t)
		for (j in 1:ncol(p)) {
			o <- !is.na(p[,j])
			p[o,j] <- p.adjust(p[o,j],method=adjust.method)
		}
		results <- new("TestResults",sign(tstat)*(p<p.value))
	},global={
		p <- as.matrix(object$p.value)
		tstat <- as.matrix(object$t)
		o <- !is.na(p)
		p[o] <- p.adjust(p[o],method=adjust.method)
		results <- new("TestResults",sign(tstat)*(p<p.value))
	},hierarchical={
		if(any(is.na(object$F.p.value))) stop("Can't handle NA p-values yet")
		sel <- p.adjust(object$F.p.value,method=adjust.method) < p.value
		i <- sum(sel,na.rm=TRUE)
		n <- sum(!is.na(sel))
		a <- switch(adjust.method,
			none=1,
			bonferroni=1/n,
			holm=1/(n-i+1),
			BH=i/n,
			BY=i/n/sum(1/(1:n))
		)
		results <- new("TestResults",array(0,dim(object$t)))
		dimnames(results) <- dimnames(object$coefficients)
		if(any(sel)) results[sel,] <- classifyTestsP(object[sel,],p.value=p.value*a,method=adjust.method)
	},nestedF={
		if(any(is.na(object$F.p.value))) stop("nestedF method can't handle NA p-values",call.=FALSE)
		sel <- p.adjust(object$F.p.value,method=adjust.method) < p.value
		i <- sum(sel,na.rm=TRUE)
		n <- sum(!is.na(sel))
		a <- switch(adjust.method,
			none=1,
			bonferroni=1/n,
			holm=1/(n-i+1),
			BH=i/n,
			BY=i/n/sum(1/(1:n))
		)
		results <- new("TestResults",array(0,dim(object$t)))
		dimnames(results) <- dimnames(object$coefficients)
		if(any(sel)) results[sel,] <- classifyTestsF(object[sel,],p.value=p.value*a)
	})
	if(lfc>0) {
		if(is.null(object$coefficients))
			warning("lfc ignored because coefficients not found")
		else
			results@.Data <- results@.Data * (abs(object$coefficients)>lfc)
	}
	results
}

classifyTestsF <- function(object,cor.matrix=NULL,df=Inf,p.value=0.01,fstat.only=FALSE) {
#	Use F-tests to classify vectors of t-test statistics into outcomes
#	Gordon Smyth
#	20 Mar 2003.  Last revised 26 June 2004.

#	Method intended for MAList objects but accept unclassed lists as well
	if(is.list(object)) {
		if(is.null(object$t)) stop("tstat cannot be extracted from object")
		if(is.null(cor.matrix) && !is.null(object$cov.coefficients)) cor.matrix <- cov2cor(object$cov.coefficients)
		if(missing(df) && !is.null(object$df.prior) && !is.null(object$df.residual)) df <- object$df.prior+object$df.residual
		tstat <- as.matrix(object$t)
	} else {
		tstat <- as.matrix(object)
	}
	ngenes <- nrow(tstat)
	ntests <- ncol(tstat)
	if(ntests == 1) {
		if(fstat.only) {
			fstat <- tstat^2
			attr(fstat,"df1") <- 1
			attr(fstat,"df2") <- df
			return(fstat)
		} else {
			p <- 2 * pt(abs(tstat), df, lower.tail=FALSE)
			return(new("TestResults", sign(tstat) * (p < p.value) ))
		}
	}

#	cor.matrix is estimated correlation matrix of the coefficients
#	and also the estimated covariance matrix of the t-statistics
	if(is.null(cor.matrix)) {
		r <- ntests
		Q <- diag(r)/sqrt(r)
	} else {
		E <- eigen(cor.matrix,symmetric=TRUE)
		r <- sum(E$values/E$values[1] > 1e-8)
		Q <- .matvec( E$vectors[,1:r], 1/sqrt(E$values[1:r]))/sqrt(r)
	}

#	Return overall moderated F-statistic only
	if(fstat.only) {
		fstat <- drop( (tstat %*% Q)^2 %*% array(1,c(r,1)) )
		attr(fstat,"df1") <- r
		attr(fstat,"df2") <- df
		return(fstat)
	}

#	Return TestResults matrix
	qF <- qf(p.value, r, df, lower.tail=FALSE)
	if(length(qF)==1) qF <- rep(qF,ngenes) 
	result <- matrix(0,ngenes,ntests,dimnames=dimnames(tstat))
	if(is.null(colnames(tstat)) && !is.null(colnames(contrasts))) colnames(result) <- colnames(contrasts)
	for (i in 1:ngenes) {
		x <- tstat[i,]
		if(any(is.na(x)))
			result[i,] <- NA
		else
			if( crossprod(crossprod(Q,x)) > qF[i] ) {
				ord <- order(abs(x),decreasing=TRUE)
				result[i,ord[1]] <- sign(x[ord[1]])
				for (j in 2:ntests) {
					bigger <- ord[1:(j-1)]
					x[bigger] <- sign(x[bigger]) * abs(x[ord[j]])
					if( crossprod(crossprod(Q,x)) > qF[i] )
						result[i,ord[j]] <- sign(x[ord[j]])
					else
						break
				}
			}
	}
	new("TestResults",result)
}

FStat <- function(object,cor.matrix=NULL)
#	Compute overall F-tests given a matrix of t-statistics
#	Gordon Smyth
#	24 February 2004.  Last modified 21 July 2004.
{
	m <- as.list(match.call())
	m[[1]] <- as.name("classifyTestsF")
	m$fstat.only <- TRUE
	eval(as.call(m))
}

classifyTestsT <- function(object,t1=4,t2=3) {
#	TestResults by rows for a matrix of t-statistics using step-down cutoffs
#	Gordon Smyth
#	1 July 2003.  Last modified 25 Feb 2004.

#	Method intended for MAList objects but accept unclassed lists as well
	if(is.list(object)) {
		if(is.null(object$t)) stop("tstat cannot be extracted from object")
		tstat <- object$t
	} else {
		tstat <- object
	}
	if(is.null(dim(tstat))) dim(tstat) <- c(1,length(tstat))
	result <- apply(tstat,1,function(x) any(abs(x)>t1,na.rm=TRUE)) * sign(tstat)*(abs(tstat)>t2)
	new("TestResults",result)
}

classifyTestsP <- function(object,df=Inf,p.value=0.05,method="holm") {
#	TestResults by rows for a matrix t-statistics using adjusted p-values
#	Gordon Smyth
#	12 July 2003.  Last modified 23 March 2004.

#	Method intended for MAList objects but accept unclassed lists as well
	if(is.list(object)) {
		if(is.null(object$t)) stop("tstat cannot be extracted from object")
		tstat <- object$t
		if(!is.null(object$df.residual)) df <- object$df.residual
		if(!is.null(object$df.prior)) df <- df+object$df.prior
	} else {
		tstat <- object
	}
	if(is.null(dim(tstat))) dim(tstat) <- c(1,length(tstat))
	ngenes <- nrow(tstat)
	P <- 2*pt(-abs(tstat),df=df)
	result <- tstat
	for (i in 1:ngenes) {
		P[i,] <- p.adjust(P[i,],method=method)
		result[i,] <- sign(tstat[i,])*(P[i,]<p.value)
	}
	new("TestResults",result)
}

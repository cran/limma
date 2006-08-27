#  DIFFERENTIAL EXPRESSION

eBayes <- function(fit,proportion=0.01,stdev.coef.lim=c(0.1,4)) {
#	Empirical Bayes statistics to select differentially expressed genes
#	Object orientated version
#	Gordon Smyth
#	4 August 2003.  Last modified 24 August 2005.

	eb <- ebayes(fit=fit,proportion=proportion,stdev.coef.lim=stdev.coef.lim)
	fit$df.prior <- eb$df.prior
	fit$s2.prior <- eb$s2.prior
	fit$var.prior <- eb$var.prior
	fit$proportion <- proportion
	fit$s2.post <- eb$s2.post
	fit$t <- eb$t
	fit$p.value <- eb$p.value
	fit$lods <- eb$lods
	if(!is.null(fit$design) && is.fullrank(fit$design)) {
		F.stat <- classifyTestsF(fit,fstat.only=TRUE)
		fit$F <- as.vector(F.stat)
		df1 <- attr(F.stat,"df1")
		df2 <- attr(F.stat,"df2")
		if(df2[1] > 1e6) # Work around bug in R 2.1
			fit$F.p.value <- pchisq(df1*fit$F,df1,lower.tail=FALSE)
		else
			fit$F.p.value <- pf(fit$F,df1,df2,lower.tail=FALSE)
	}
	fit
}

ebayes <- function(fit,proportion=0.01,stdev.coef.lim=c(0.1,4)) {
#	Empirical Bayes statistics to select differentially expressed genes
#	Gordon Smyth
#	8 Sept 2002.  Last revised 29 May 2004.

	coefficients <- fit$coefficients
	stdev.unscaled <- fit$stdev.unscaled
	sigma <- fit$sigma
	df.residual <- fit$df.residual
	if(is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || is.null(df.residual)) stop("No data, or argument is not a valid lmFit object")
	if(all(df.residual==0)) stop("No residual degrees of freedom in linear model fits")
	if(all(!is.finite(sigma))) stop("No finite residual standard deviations")

#	Moderated t-statistic
	out <- squeezeVar(sigma^2, df.residual)
	out$s2.prior <- out$var.prior
	out$s2.post <- out$var.post
	out$var.prior <- out$var.post <- NULL
	df.total <- df.residual + out$df.prior
	out$t <- coefficients / stdev.unscaled / sqrt(out$s2.post)
	out$p.value <- 2*pt(-abs(out$t),df=df.total)

#	B-statistic
	var.prior.lim <- stdev.coef.lim^2/out$s2.prior
	out$var.prior <- tmixture.matrix(out$t,stdev.unscaled,df.total,proportion,var.prior.lim)
	if(any(is.na(out$var.prior))) {
		out$var.prior[ is.na(out$var.prior) ] <- 1/out$s2.prior
		warning("Estimation of var.prior failed - set to default value")
	}
	r <- rep(1,NROW(out$t)) %o% out$var.prior
	r <- (stdev.unscaled^2+r) / stdev.unscaled^2
	t2 <- out$t^2
	if(out$df.prior > 10^6)
		kernel <- t2*(1-1/r)/2
	else
		kernel <- (1+df.total)/2*log((t2+df.total) / (t2/r+df.total))
	out$lods <- drop( log(proportion/(1-proportion))-log(r)/2+kernel )
	out
}

tmixture.matrix <- function(tstat,stdev.unscaled,df,proportion,v0.lim=NULL) {
#	Estimate the prior variance of the coefficients for DE genes
#	Gordon Smyth
#	18 Nov 2002. Last modified 12 Dec 2003.

	tstat <- as.matrix(tstat)
	stdev.unscaled <- as.matrix(stdev.unscaled)
	if(any(dim(tstat) != dim(stdev.unscaled))) stop("Dims of tstat and stdev.unscaled don't match")
	if(!is.null(v0.lim)) if(length(v0.lim) != 2) stop("v0.lim must have length 2")
	ncoef <- ncol(tstat)
	v0 <- rep(0,ncoef)
	for (j in 1:ncoef) v0[j] <- tmixture.vector(tstat[,j],stdev.unscaled[,j],df,proportion,v0.lim)	
	v0
}

tmixture.vector <- function(tstat,stdev.unscaled,df,proportion,v0.lim=NULL) {
#	Estimate scale factor in mixture of two t-distributions
#	tstat is assumed to follow (v0+v1)/v1*t(df) with probability proportion and t(df) otherwise
#	v1 is stdev.unscaled^2 and v0 is to be estimated
#	Gordon Smyth
#	18 Nov 2002.  Last modified 13 Dec 2003.

	if(any(is.na(tstat))) {
		o <- !is.na(tstat)
		tstat <- tstat[o]
		stdev.unscaled <- stdev.unscaled[o]
		df <- df[o]
	}
	ngenes <- length(tstat)
	ntarget <- ceiling(proportion/2*ngenes)
	if(ntarget < 1) return(NA)

#	If ntarget is v small, ensure p at least matches selected proportion
#	This ensures ptarget < 1
	p <- max(ntarget/ngenes,proportion)

	tstat <- abs(tstat)
	ttarget <- quantile(tstat,(ngenes-ntarget)/(ngenes-1))
	top <- (tstat >= ttarget)
	tstat <- tstat[top]
	v1 <- stdev.unscaled[top]^2
	df <- df[top]
	r <- ntarget-rank(tstat)+1
	p0 <- pt(-tstat,df=df)
	ptarget <- ( (r-0.5)/2/ngenes - (1-p)*p0 ) / p
	pos <- ptarget > p0
	v0 <- rep(0,ntarget)
	if(any(pos)) {
		qtarget <- qt(ptarget[pos],df=df[pos])
		v0[pos] <- v1[pos]*((tstat[pos]/qtarget)^2-1)
	}
	if(!is.null(v0.lim)) v0 <- pmin(pmax(v0,v0.lim[1]),v0.lim[2])
	mean(v0)
}

fitFDist <- function(x,df1) {
#	Moment estimation of the parameters of a scaled F-distribution
#	The first degrees of freedom is given
#	Gordon Smyth
#	8 Sept 2002.  Last revised 6 April 2006.

#	Remove missing or infinite values and zero degrees of freedom
	o <- is.finite(x) & is.finite(df1) & (x >= 0) & (df1 > 0)
	if(any(!o)) {
		x <- x[o]
		df1 <- df1[o]
	}
	n <- length(x)
	if(n==0) return(list(scale=NA,df2=NA))

#	Avoid exactly zero values
	m <- median(x)
	if(m==0) {
		warning("More than half of residual variances are exactly zero: eBayes unreliable")
		m <- 1
	} else {
		if(any(x==0)) warning("Zero sample variances detected, have been offset",call.=FALSE)
	}
	x <- pmax(x, 1e-5 * median(x))

#	Better to work on with log(F)
	z <- log(x)
	e <- z-digamma(df1/2)+log(df1/2)
	emean <- mean(e)
	evar <- mean(n/(n-1)*(e-emean)^2-trigamma(df1/2))
	if(evar > 0) {
		df2 <- 2*trigammaInverse(evar)
		s20 <- exp(emean+digamma(df2/2)-log(df2/2))
	} else {
		df2 <- Inf
		s20 <- exp(emean)
	}
	list(scale=s20,df2=df2)
}

trigammaInverse <- function(x) {
#	Solve trigamma(y) = x for y
#	Gordon Smyth
#	8 Sept 2002.  Last revised 12 March 2004.

#	Non-numeric or zero length input
	if(!is.numeric(x)) stop("Non-numeric argument to mathematical function")
	if(length(x)==0) return(numeric(0))

#	Treat out-of-range values as special cases
	omit <- is.na(x)
	if(any(omit)) {
		y <- x
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 0)
	if(any(omit)) {
		y <- x
		y[omit] <- NaN
		warning("NaNs produced")
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x > 1e7)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/sqrt(x[omit])
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}
	omit <- (x < 1e-6)
	if(any(omit)) {
		y <- x
		y[omit] <- 1/x[omit]
		if(any(!omit)) y[!omit] <- Recall(x[!omit])
		return(y)
	}

#	Newton's method
#	1/trigamma(y) is convex, nearly linear and strictly > y-0.5,
#	so iteration to solve 1/x = 1/trigamma is monotonically convergent
	y <- 0.5+1/x
	iter <- 0
	repeat {
		iter <- iter+1
		tri <- trigamma(y)
		dif <- tri*(1-tri/x)/psigamma(y,deriv=2)
		y <- y+dif
		if(max(-dif/y) < 1e-8) break
		if(iter > 50) {
			warning("Iteration limit exceeded")
			break
		}
	}
	y
}

qqt <- function(y,df=Inf,ylim=range(y),main="Student's t Q-Q Plot",xlab="Theoretical Quantiles",ylab="Sample Quantiles",plot.it=TRUE,...)
{
#	Student's t probability plot
#	Gordon Smyth
#	3 Oct 2002

    y <- y[!is.na(y)]
    if(0 == (n <- length(y))) stop("y is empty")
    x <- qt(ppoints(n),df=df)[order(order(y))]
    if (plot.it) plot(x,y,main=main,xlab=xlab,ylab=ylab,ylim=ylim,...)
    invisible(list(x=x,y=y))
}

squeezeVar <- function(var, df)
#	Empirical Bayes posterior variances
#	Gordon Smyth
#	2 March 2004
{
	n <- length(var)
	if(n == 0) stop("var is empty")
	if(n == 1) return(list(var.post=var,var.prior=var,df.prior=0))
	if(length(df)==1) { 
		df <- rep.int(df,n)
	} else {
		if(length(df) != n) stop("lengths differ")
	}
	out <- fitFDist(var, df1=df)
	if(is.null(out$df2) || is.na(out$df2)) stop("Could not estimate prior df")
	out$var.prior <- out$scale
	out$df.prior <- out$df2
	out$df2 <- out$scale <- NULL
	df.total <- df + out$df.prior
	if(out$df.prior == Inf)
		out$var.post <- rep.int(out$var.prior,n)
	else {
		var[df==0] <- 0 # guard against missing or infinite values
		out$var.post <- (df*var + out$df.prior*out$var.prior) / df.total
	}
	out
}


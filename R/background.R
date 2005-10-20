#  BACKGROUND.R

#  BACKGROUND CORRECTION

backgroundCorrect <- function(RG, method="subtract", offset=0, printer=RG$printer, verbose=TRUE) {
#	Apply background correction to microarray data
#	Gordon Smyth
#	12 April 2003.  Last modified 8 Oct 2004.

	if(is.null(RG$Rb) != is.null(RG$Gb)) stop("Background values exist for one channel but not the other")
	method <- match.arg(method, c("none","subtract","half","minimum","movingmin","edwards","normexp"))
	if(is.null(RG$Rb) && is.null(RG$Gb)) method <- "none"
	switch(method,
	subtract={
		RG$R <- RG$R-RG$Rb
		RG$G <- RG$G-RG$Gb
	},
	half={
		RG$R <- pmax(RG$R-RG$Rb, 0.5)
		RG$G <- pmax(RG$G-RG$Gb, 0.5)
	},
	minimum={
		RG$R <- as.matrix(RG$R - RG$Rb)
		RG$G <- as.matrix(RG$G - RG$Gb)
		for (slide in 1:ncol(RG$R)) {
			i <- RG$R[,slide] < 1e-18
			if(any(i,na.rm=TRUE)) {
				m <- min(RG$R[!i,slide],na.rm=TRUE)
				RG$R[i,slide] <- m/2
			}
			i <- RG$G[,slide] < 1e-18
			if(any(i,na.rm=TRUE)) {
				m <- min(RG$G[!i,slide],na.rm=TRUE)
				RG$G[i,slide] <- m/2
			}
		}
	},
	movingmin={
		RG$R <- RG$R-ma3x3.spottedarray(RG$Rb,printer=printer,FUN=min,na.rm=TRUE)
		RG$G <- RG$G-ma3x3.spottedarray(RG$Gb,printer=printer,FUN=min,na.rm=TRUE)
	},
	edwards={
#		Log-linear interpolation for dull spots as in Edwards (2003).
#		The threshold values (delta) are chosen such that the number of
#		spots with (0 < R-Rb < delta) is f=10% of the number spots
#		with (R-Rb <= 0) for each channel and array.
#		Note slight change from Edwards (2003).
		one <- matrix(1,NROW(RG$R),1)
		delta.vec <- function(d, f=0.1) {
			quantile(d, mean(d<1e-16,na.rm=TRUE)*(1+f), na.rm=TRUE)
		}
		sub <- as.matrix(RG$R-RG$Rb)
		delta <- one %*% apply(sub, 2, delta.vec)
		RG$R <- ifelse(sub < delta, delta*exp(1-(RG$Rb+delta)/RG$R), sub)
		sub <- as.matrix(RG$G-RG$Gb)
		delta <- one %*% apply(sub, 2, delta.vec)
		RG$G <- ifelse(sub < delta, delta*exp(1-(RG$Gb+delta)/RG$G), sub)
	},
	normexp={
	for (j in 1:ncol(RG$R)) {
		out <- normexp.fit(foreground=RG$G[,j],background=RG$Gb[,j])
		RG$G[,j] <- normexp.signal(mu=out$beta+RG$Gb[,j],out$sigma,out$alpha,foreground=RG$G[,j])
		out <- normexp.fit(foreground=RG$R[,j],background=RG$Rb[,j])
		RG$R[,j] <- normexp.signal(mu=out$beta+RG$Rb[,j],out$sigma,out$alpha,foreground=RG$R[,j])
		if(verbose) cat("Corrected array",j,"\n")
	}
	})
	RG$Rb <- NULL
	RG$Gb <- NULL
	if(offset) {
		RG$R <- RG$R+offset
		RG$G <- RG$G+offset
	}
	new("RGList",unclass(RG))
}

ma3x3.matrix <- function(x,FUN=mean,na.rm=TRUE,...)
#	2-dimensional moving average for 3x3 blocks
#	Gordon Smyth
#	11 April 2004
{
#	Pad out x with NA so that original values have 8 neighbors
	d1 <- nrow(x)
	d2 <- ncol(x)
	y <- matrix(NA,d1+2,d2+2)
	y[1+(1:d1),1+(1:d2)] <- x

#	Index vector for original values
	i <- 1:length(y)
	dim(i) <- dim(y)
	i <- i[1+(1:d1),1+(1:d2)]
	dim(i) <- NULL

#	Rows are original obs, columns are neighbors
	x <- matrix(x,d1*d2,9)
	ry <- nrow(y)
	x[,1] <- y[i-ry-1]
	x[,2] <- y[i-ry]
	x[,3] <- y[i-ry+1]
	x[,4] <- y[i-1]
	x[,6] <- y[i+1]
	x[,7] <- y[i+ry-1]
	x[,8] <- y[i+ry]
	x[,9] <- y[i+ry+1]

	y <- apply(x,MARGIN=1,FUN=FUN,na.rm=na.rm,...)
	dim(y) <- c(d1,d2)
	y
}

ma3x3.spottedarray <- function(x,printer,FUN=mean,na.rm=TRUE,...)
#	Gordon Smyth
#	11 April 2004
{
	x <- as.matrix(x)
	narrays <- ncol(x)
	gr <- printer$ngrid.r
	gc <- printer$ngrid.c
	sr <- printer$nspot.r
	sc <- printer$nspot.c
	dim(x) <- c(sc, sr, gc, gr, narrays)
	x <- aperm(x, perm = c(2, 4, 1, 3, 5))
	dim(x) <- c(gr * sr, gc * sc, narrays)
	for (j in 1:narrays) x[,,j] <- ma3x3.matrix(x[,,j],FUN=FUN,na.rm=TRUE,...)
	dim(x) <- c(sr, gr, sc, gc, narrays)
	x <- aperm(x, perm = c(3, 1, 4, 2, 5))
	dim(x) <- c(sc*sr*gc*gr, narrays)
	x
}

#  NORMAL + EXPONENTIAL ADAPTIVE MODEL

normexp.signal <- function(mu,sigma,alpha,foreground)
#	Expected value of signal given foreground in normal + exponential model
#	Gordon Smyth
#	24 Aug 2002. Last modified 26 Sept 2005.
{
	if(alpha <= 0) stop("alpha must be positive")
	if(sigma <= 0) stop("sigma must be positive")
	mu.sf <- foreground-mu-sigma^2/alpha
	signal <- mu.sf + sigma^2 * exp(dnorm(0,mean=mu.sf,sd=sigma,log=TRUE) - pnorm(0,mean=mu.sf,sd=sigma,lower.tail=FALSE,log=TRUE))
	o <- !is.na(signal)
	if(any(signal[o]<0)) {
		warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensity to small value")
		signal[o] <- pmax(signal[o],1e-6)
	}
	signal
}

normexp.fit <- function(foreground,background=0,trace=0)
#	Fit background=normal + signal=exponential model using BFGS.
#	Gordon Smyth and Jeremy Silver
#	24 Aug 2002. Last modified 6 October 2005.
{
	x <- foreground-background
	isna <- is.na(x)
	if(any(isna)) x <- x[!isna]
	if(length(x)<4) stop("Not enough data: need at least 4 non-missing corrected intensities")
	
#	Starting values for parameters beta, alpha and sigma
	q <- quantile(x, c(0,0.05,0.1,1), na.rm = TRUE, names = FALSE)
	if(q[1]==q[4]) return(list(beta=q[1],sigma=1,alpha=1,m2loglik=NA,convergence=0))
	if(q[2] > q[1]) {
		beta <- q[2]
	} else {
		if(q[3] > q[1]) {
			beta <- q[3]
		} else {
			beta <- q[1] + 0.05*(q[4]-q[1])
		}
	}
	sigma <- sqrt(mean((x[x<beta]-beta)^2, na.rm = TRUE))
	alpha <- mean(x,na.rm = TRUE) - beta
	if(alpha <= 0) alpha <- 1e-6

	Results <- optim(par=c(beta,log(sigma),log(alpha)), fn=normexp.m2loglik, gr=normexp.grad, method=c("BFGS"), control=list(trace=trace), foreground=x, background=0)
	list(beta=Results$par[1], sigma=exp(Results$par[2]), alpha=exp(Results$par[3]),m2loglik=Results$value, convergence=Results$convergence)
}

normexp.grad <- function(theta,foreground,background=0)
#	Gradient of norm-exp log-likelihood (summed over all spots)
#	Jeremy Silver.
#	21 Jan 2005. Last modified 18 June 2005.
{
	beta<-theta[1]
	logsigma<-theta[2]
	logalpha<-theta[3]
	mu <- beta + background
	mu.sf <- foreground - mu - exp(2*logsigma - logalpha)

	dlogdbeta <- -2 * sum(exp(-logalpha) - exp(dnorm(0,mu.sf,exp(logsigma),log = TRUE) - pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE,log.p = TRUE)))
	dlogdlogsigma <- -2 * sum(exp(2*logsigma - 2*logalpha) - (2*exp(2*logsigma - logalpha) + mu.sf)*exp( dnorm(0,mu.sf,exp(logsigma),log = TRUE) - pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE,log.p = TRUE)))
	dlogdlogalpha <- -2 * sum(-1 +(foreground - mu)*exp( -logalpha) - exp(2*logsigma - 2*logalpha) + exp(2*logsigma - logalpha+dnorm(0,mu.sf,exp(logsigma),log = TRUE) - pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE,log.p = TRUE)))

	c(dlogdbeta,dlogdlogsigma,dlogdlogalpha)
}

normexp.m2loglik <- function(theta,foreground,background=0)
#	Minus twice the norm-exp log-likelihood (summed over all spots)
#	Jeremy Silver and Gordon Smyth
#	24 Aug 2002. Last modified 6 October 2005.
{
	beta<-theta[1]
	logsigma<-theta[2]
	logalpha<-theta[3]
	mu <- beta + background
	mu.sf <- foreground - mu - exp(2*logsigma - logalpha)
	
	m2loglik <- -2*sum(-logalpha - (foreground - mu)/exp(logalpha) + 0.5*exp(2*logsigma - 2*logalpha) + pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE, log.p=TRUE))
	max(min(m2loglik,.Machine$double.xmax),.Machine$double.xmin)
}


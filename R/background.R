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
		out <- fit.normexp(foreground=RG$G[,j],background=RG$Gb[,j])
		RG$G[,j] <- signal.normexp(mu=out$beta+RG$Gb[,j],out$sigma,out$alpha,foreground=RG$G[,j])
		out <- fit.normexp(foreground=RG$R[,j],background=RG$Rb[,j])
		RG$R[,j] <- signal.normexp(mu=out$beta+RG$Rb[,j],out$sigma,out$alpha,foreground=RG$R[,j])
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

fit.normexp <- function(foreground,background=NULL,background.matrix=NULL,trace=0,beta.start=NULL) {
#	Fit background=normal + signal=exponential model.
#	Gordon Smyth
#	24 Aug 2002.  Last modified 19 Nov 2004.

	f <- foreground

#	Starting values
	mu <- quantile(f,0.03,na.rm=TRUE,names=FALSE)
	if(is.null(background) && is.null(background.matrix)) {
		nbeta <- 1
		beta <- mu
	} else {
		if(!is.null(background)) {
			nbeta <- 1
			beta <- mu - mean(background,na.rm=TRUE)
		} else {
			nbeta <- ncol(background.matrix)
			beta <- rep(0,nbeta)
			beta[1] <- mu
			if(!is.null(beta.start)) beta <- beta.start
		}
	}
	sigma <- sqrt(mean((f[f<mu]-mu)^2,na.rm=TRUE))
	if(!is.infinite(sigma) || sigma < 1) sigma <- 1
	alpha <- mean(f,na.rm=TRUE) - mu
	theta <- c(beta,log(sigma),log(alpha))
	
#	Nelder-Mead optimization
	out <- optim(theta,m2loglik.normexp,control=list(trace=trace),foreground=f,background=background,background.matrix=background.matrix)

	if(out$convergence==1) warning("optim iteration limit reached")
	if(out$convergence==10) stop("optim failure: degeneracy of the Nelder-Mead simplex")
	list(beta=out$par[1:nbeta],sigma=exp(out$par[nbeta+1]),alpha=exp(out$par[nbeta+2]),m2loglik=out$value,convergence=out$convergence)
}

m2loglik.normexp <- function(theta,foreground,background=NULL,background.matrix=NULL) {
#	Marginal log-likelihood of foreground values for normal + exponential model.
#	Gordon Smyth
#	24 Aug 2002.

	if(is.null(background) && is.null(background.matrix)) {
		nbeta <- 1
		mu <- theta[1]
	} else {
		if(!is.null(background)) {
			nbeta <- 1
			mu <- theta[1]+background
		} else {
			nbeta <- ncol(background.matrix)
			mu <- background.matrix %*% theta[1:nbeta]
		}
	}
	sigma <- exp(theta[nbeta+1])
	alpha <- exp(theta[nbeta+2])
	mu.sf <- foreground-mu-sigma^2/alpha
	-2*sum(-log(alpha) - (foreground-mu)/alpha + sigma^2/alpha^2/2 + pnorm(0,mean=mu.sf,sd=sigma,lower.tail=FALSE,log.p=TRUE))
}

signal.normexp <- function(mu,sigma,alpha,foreground) {
#	Expected value of signal given foreground in normal + exponential model
#	Gordon Smyth
#	24 Aug 2002.

	mu.sf <- foreground-mu-sigma^2/alpha
	mu.sf + sigma^2 * dnorm(0,mean=mu.sf,sd=sigma) / pnorm(0,mean=mu.sf,sd=sigma,lower.tail=FALSE)
}

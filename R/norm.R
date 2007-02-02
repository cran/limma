#	LOESS FUNCTIONS

#loessFitnew <- function(y, x, weights=NULL, span=0.3, bin=0.01/(2-is.null(weights)), iterations=4) {
##	Fast loess fit for simple x and y
##	This function uses lowess if no weights and loess otherwise.
##	It is intended to give a streamlined common interface to the two functions.
##	Gordon Smyth
##	28 June 2003.  Last revised 27 Feb 2005.
#
#	n <- length(y)
#	out <- list(fitted=rep(NA,n),residuals=rep(NA,n))
#	obs <- is.finite(y) & is.finite(x)
#	xobs <- x[obs]
#	yobs <- y[obs]
#	nobs <- length(yobs)
#	if(nobs==0) return(out)
#	if(is.null(weights)) {
#		o <- order(xobs)
#		oo <- order(o)
#		iter <- iterations-1
#		delta = bin * diff(range(xobs)) 
#		smoothy <- .C("lowess", x = as.double(xobs[o]), as.double(yobs[o]), 
#			nobs, as.double(span), as.integer(iter), as.double(delta), 
#			y = double(nobs), double(nobs), double(nobs), PACKAGE = "base")$y[oo]
#		out$fitted[obs] <- smoothy
#		out$residuals[obs] <- yobs-smoothy
#	} else {
#		weights[is.na(weights) | weights<0] <- 0
#		wt0 <- obs & weights==0
#		if(any(wt0)) {
#			obs[wt0] <- FALSE
#			xobs <- x[obs]
#			yobs <- y[obs]
#			nobs <- length(yobs)
#			if(nobs==0) return(out)
#		}
#		wobs <- weights[obs]
#		if(nobs < 4+1/span) {
#			fit <- lm.wfit(cbind(1,xobs),yobs,wobs)
#		} else {
##			Suppress warning "k-d tree limited by memory"
#			oldopt <- options(warn=-1)
#			on.exit(options(oldopt))
# This is now quite fast
#			fit <- loess(yobs~xobs,weights=wobs,span=span,degree=1,family="symmetric",cell=bin/span,iterations=iterations,trace.hat="approximate")
#		}
#		out$fitted[obs] <- fitted(fit)
#		out$residuals[obs] <- residuals(fit)
#		if(any(wt0)) {
#			out$fitted[wt0] <- predict(fit,newdata=x[wt0])
#			out$residuals[wt0] <- y[wt0]-out$fitted(
#		}
#	}
#	out
#}

loessFit <- function(y, x, weights=NULL, span=0.3, bin=0.01/(2-is.null(weights)), iterations=4) {
#	Fast loess fit for simple x and y
#	This function uses lowess if no weights and loess otherwise.
#	It is intended to give a streamlined common interface to the two functions.
#	Gordon Smyth
#	28 June 2003.  Last revised 16 Feb 2004.

	n <- length(y)
	out <- list(fitted=rep(NA,n),residuals=rep(NA,n))
	obs <- is.finite(y) & is.finite(x)
	xobs <- x[obs]
	yobs <- y[obs]
	nobs <- length(yobs)
	if(nobs==0) return(out)
	if(is.null(weights)) {
		o <- order(xobs)
		oo <- order(o)
		iter <- iterations-1
		delta = bin * diff(range(xobs)) 
		smoothy <- .C("lowess", x = as.double(xobs[o]), as.double(yobs[o]), 
			nobs, as.double(span), as.integer(iter), as.double(delta), 
			y = double(nobs), double(nobs), double(nobs), PACKAGE = "base")$y[oo]
		out$fitted[obs] <- smoothy
		out$residuals[obs] <- yobs-smoothy
	} else {
		weights[!is.finite(weights)] <- 0
		wobs <- weights[obs]
		if(sum(wobs>0) < 4+1/span) {
			fit <- lm.wfit(cbind(1,xobs),yobs,wobs)
		} else {
#			Suppress warning "k-d tree limited by memory"
#			oldopt <- options(warning.expression=expression())
			oldopt <- options(warn=-1)
			on.exit(options(oldopt))
			fit <- .vsimpleLoess(y=yobs, x=xobs, weights=wobs, span=span, degree=1,
				cell=bin/span, iterations=iterations)
		}
		out$fitted[obs] <- fit$fitted
		out$residuals[obs] <- fit$residuals
	}
	out
}

.vsimpleLoess <- function (y, x, weights, span=0.75, degree=2, cell=0.2, iterations=1)
#	Cut-down version of simpleLoess from modreg package
#	Gordon Smyth
#	28 June 2003.  Last modified 7 March 2004.
{
	statistics <- "none"
	surface <- "interpolate"
	surf.stat <- paste(surface, statistics, sep = "/")
	D <- 1
	N <- NROW(x)
	if (!N || !D) stop("invalid `x'")
	if (!length(y)) stop("invalid `y'")
	x <- as.matrix(x)
	max.kd <- max(N, 200)
	robust <- rep(1, N)
	sum.drop.sqr <- 0
	sum.parametric <- 0
	nonparametric <- 1
	order.parametric <- 1
	order.drop.sqr <- 2
	if (iterations) 
		for (j in 1:iterations) {
			robust <- weights * robust
			z <- .C("loess_raw", as.double(y), as.double(x), 
				as.double(weights), as.double(robust), as.integer(D), 
				as.integer(N), as.double(span), as.integer(degree), 
				as.integer(nonparametric), as.integer(order.drop.sqr), 
				as.integer(sum.drop.sqr), as.double(span * cell), 
				as.character(surf.stat), fitted.values = double(N), 
				parameter = integer(7), a = integer(max.kd), 
				xi = double(max.kd), vert = double(2 * D), vval = double((D + 
				1) * max.kd), diagonal = double(N), trL = double(1), 
				delta1 = double(1), delta2 = double(1), as.integer(surf.stat == 
				"interpolate/exact"), PACKAGE = "stats")
			fitted.residuals <- y - z$fitted.values
			if (j < iterations) 
				robust <- .Fortran("lowesw", as.double(fitted.residuals), 
				as.integer(N), robust = double(N), double(N), 
				PACKAGE = "stats")$robust
		}
	list(fitted = z$fitted.values, residuals = fitted.residuals)
}

#  WITHIN ARRAY NORMALIZATION

MA.RG <- function(object, bc.method="subtract", offset=0) {
#	Convert RGList to MAList
#	Gordon Smyth
#	2 March 2003.  Last revised 9 Dec 2004.

	if(is.null(object$R) || is.null(object$G)) stop("Object doesn't contain R and G components")
	object <- backgroundCorrect(object, method=bc.method, offset=offset)
	R <- object$R
	G <- object$G

#	Log
	R[R <= 0] <- NA
	G[G <= 0] <- NA
	R <- log(R,2)
	G <- log(G,2)
	
#	Minus and Add
	object$R <- object$G <- object$Rb <- object$Gb <- object$other <- NULL
	object$M <- as.matrix(R-G)
	object$A <- as.matrix((R+G)/2)
	new("MAList",unclass(object))
}

RG.MA <- function(object) {
#	Convert MAList to RGList
#	Gordon Smyth
#	5 September 2003.  Last modified 9 Dec 2004.

	object$R <- 2^( object$A+object$M/2 )
	object$G <- 2^( object$A-object$M/2 )
	object$M <- NULL
	object$A <- NULL
	new("RGList",unclass(object))
}

normalizeWithinArrays <- function(object,layout=object$printer,method="printtiploess",weights=object$weights,span=0.3,iterations=4,controlspots=NULL,df=5,robust="M",bc.method="subtract",offset=0)
#	Within array normalization
#	Gordon Smyth
#	2 March 2003.  Last revised 2 September 2006.
{
#	Check input arguments
#	and get non-intensity dependent methods out of the way
	if(!is(object,"MAList")) object <- MA.RG(object,bc.method=bc.method,offset=offset)
	choices <- c("none","median","loess","printtiploess","composite","control","robustspline")
	method <- match.arg(method,choices)
	if(method=="none") return(object)
	if(is.vector(object$M)) object$M <- as.matrix(object$M)
	nprobes <- nrow(object$M)
	narrays <- ncol(object$M)
	if(!is.null(weights)) weights <- asMatrixWeights(weights,dim=c(nprobes,narrays))
	if(method=="median") {
		if(is.null(weights))
			for (j in 1:narrays) object$M[,j] <- object$M[,j] - median(object$M[,j],na.rm=TRUE)
		else
			for (j in 1:narrays) object$M[,j] <- object$M[,j] - weighted.median(object$M[,j],weights[,j],na.rm=TRUE)
		return(object)
	}

#	All remaining methods use regression of M-values on A-values
	if(is.vector(object$A)) object$A <- as.matrix(object$A)
	if(nrow(object$A) != nprobes) stop("row dimension of A doesn't match M")
	if(ncol(object$A) != narrays) stop("col dimension of A doesn't match M")
	switch(method,
		loess = {
			for (j in 1:narrays) {
				y <- object$M[,j]
				x <- object$A[,j]
				w <- weights[,j]
				object$M[,j] <- loessFit(y,x,w,span=span,iterations=iterations)$residuals
			}
		},
		printtiploess = {
			if(is.null(layout)) stop("Layout argument not specified")
			ngr <- layout$ngrid.r
			ngc <- layout$ngrid.c
			nspots <- layout$nspot.r * layout$nspot.c
			nprobes2 <- ngr*ngc*nspots
			if(nprobes2 != nprobes) stop("printer layout information does not match M row dimension")
			for (j in 1:narrays) {
				spots <- 1:nspots
				for (gridr in 1:ngr)
				for (gridc in 1:ngc) {
					y <- object$M[spots,j]
					x <- object$A[spots,j]
					w <- weights[spots,j]
					object$M[spots,j] <- loessFit(y,x,w,span=span,iterations=iterations)$residuals
					spots <- spots + nspots
				}
			}
		},
		composite = {
			if(is.null(layout)) stop("Layout argument not specified")
			if(is.null(controlspots)) stop("controlspots argument not specified")
			ntips <- layout$ngrid.r * layout$ngrid.c
			nspots <- layout$nspot.r * layout$nspot.c
			for (j in 1:narrays) {
				y <- object$M[,j]
				x <- object$A[,j]
				w <- weights[,j]
				f <- is.finite(y) & is.finite(x)
				if(!is.null(w)) f <- f & is.finite(w)
				y[!f] <- NA
				fit <- loess(y~x,weights=w,span=span,subset=controlspots,na.action=na.exclude,degree=0,surface="direct",family="symmetric",trace.hat="approximate",iterations=iterations)
				alpha <- global <- y
				global[f] <- predict(fit,newdata=x[f])
				alpha[f] <- (rank(x[f])-1) / sum(f)
				spots <- 1:nspots
				for (tip in 1:ntips) {
					y <- object$M[spots,j]
					x <- object$A[spots,j]
					w <- weights[spots,j]
					local <- loessFit(y,x,w,span=span,iterations=iterations)$fitted
					object$M[spots,j] <- object$M[spots,j] - alpha[spots]*global[spots]-(1-alpha[spots])*local
					spots <- spots + nspots
				}
			}
		},
		control = {
			if(is.null(layout)) stop("Layout argument not specified")
			if(is.null(controlspots)) stop("controlspots argument not specified")
			ntips <- layout$ngrid.r * layout$ngrid.c
			nspots <- layout$nspot.r * layout$nspot.c
			for (j in 1:narrays) {
				y <- object$M[,j]
				x <- object$A[,j]
				w <- weights[,j]
				f <- is.finite(y) & is.finite(x)
				if(!is.null(w)) f <- f & is.finite(w)
				y[!f] <- NA
				fit <- loess(y~x,weights=w,span=span,subset=controlspots,na.action=na.exclude,degree=1,surface="direct",family="symmetric",trace.hat="approximate",iterations=iterations)
				y[f] <- y[f]-predict(fit,newdata=x[f])
				object$M[,j] <- y
			}
		},
		robustspline = {
			if(is.null(layout)) stop("Layout argument not specified")
			for (j in 1:narrays)
				object$M[,j] <- normalizeRobustSpline(object$M[,j],object$A[,j],layout,df=df,method=robust)
		}
	)
	object
}

normalizeRobustSpline <- function(M,A,layout,df=5,method="M") {
#	Robust spline normalization
#	Gordon Smyth
#	27 April 2003.  Last revised 2 Feb 2007.

	require(MASS)
	require(splines)
	ngrids <- round(layout$ngrid.r * layout$ngrid.c)
	nspots <- round(layout$nspot.r * layout$nspot.c)
	if(ngrids < 1) stop("layout incorrectly specified")
	if(nspots < df+1) stop("too few spots in each print-tip block")

#	Global splines
	O <- is.finite(M) & is.finite(A)
	X <- matrix(NA,ngrids*nspots,df)
	X[O,] <- ns(A[O],df=df,intercept=TRUE)
	x <- X[O,,drop=FALSE]
	y <- M[O]
	w <- rep(1,length(y))
	s0 <- summary(rlm(x,y,weights=w,method=method),method="XtWX",correlation=FALSE)
	beta0 <- s0$coefficients[,1]
	covbeta0 <- s0$cov * s0$stddev^2

#	Case with only one tip-group
	if(ngrids <= 1) {
		M[O] <- s0$residuals
		M[!O] <- NA
		return(M)
	}

#	Tip-wise splines
	beta <- array(1,c(ngrids,1)) %*% array(beta0,c(1,df))
	covbeta <- array(0,c(ngrids,df,df))
	spots <- 1:nspots
	for (i in 1:ngrids) {
		o <- O[spots]
		y <- M[spots][o]
		if(length(y)) {
			x <- X[spots,][o,,drop=FALSE]
			r <- qr(x)$rank
			if(r<df) x <- x[,1:r,drop=FALSE]
			w <- rep(1,length(y))
			s <- summary(rlm(x,y,weights=w,method=method),method="XtWX",correlation=FALSE)
			beta[i,1:r] <- s$coefficients[,1]
			covbeta[i,1:r,1:r] <- s$cov * s$stddev^2
		}
		spots <- spots + nspots
	}

#	Empirical Bayes estimates
	res.cov <- cov(beta) - apply(covbeta,c(2,3),mean)
	a <- max(0, sum(diag(res.cov)) / sum(diag(covbeta0)) )
	Sigma0 <- covbeta0 * a 
#	Sigma0 <- covbeta0 * max(0,mean(eigen(solve(covbeta0,res.cov))$values))

#	Shrunk splines
	if(a==0) {
		M[O] <- s0$residuals
		M[!O] <- NA
	} else {
		spots <- 1:nspots
		for (i in 1:ngrids) {
			beta[i,] <- beta0 + Sigma0 %*% solve(Sigma0+covbeta[i,,],beta[i,]-beta0)
			o <- O[spots]
			x <- X[spots,][o,]
			M[spots][o] <- M[spots][o] - x %*% beta[i,]
			M[spots][!o] <- NA
			spots <- spots + nspots
		}
	}
	M
}


#  PRINTORDER

normalizeForPrintorder <- function(object,layout,start="topleft",method="loess",separate.channels=FALSE,span=0.1,plate.size=32) {
#	Pre-normalize the foreground intensities for print order
#	Gordon Smyth
#	11 Mar 2002.  Last revised 18 June 2003.

	if(is.null(object$R) || is.null(object$G)) stop("List must contain components R and G")
	start <- match.arg(start,c("topleft","topright"))
	method <- match.arg(method,c("loess","plate"))
	po <- printorder(layout,start=start)$printorder
	nslides <- NCOL(object$R)
	for (i in 1:nslides) {
		RG <- normalizeForPrintorder.rg(R=object$R[,i],G=object$G[,i],printorder=po,method=method,separate.channels=separate.channels,span=span,plate.size=plate.size)
		object$R[,i] <- RG$R
		object$G[,i] <- RG$G
	}
	object
}

normalizeForPrintorder.rg <- function(R,G,printorder,method="loess",separate.channels=FALSE,span=0.1,plate.size=32,plot=FALSE) {
#	Pre-normalize the foreground intensities for print order, given R and G for a single array.
#	Gordon Smyth
#	8 Mar 2002.  Last revised 3 January 2007.

	if(plot) ord <- order(printorder)
	Rf <- log(R,2)
	Gf <- log(G,2)
	Rf[is.infinite(Rf)] <- NA
	Gf[is.infinite(Gf)] <- NA
	if(!separate.channels) Af <- (Rf+Gf)/2
	method <- match.arg(method,c("loess","plate"))
	if(method=="plate") {
		# Correct for plate pack (usually four 384-well plates)
		plate <- 1 + (printorder-0.5) %/% plate.size
		hubermu <- function(...) huber(...)$mu
		if(separate.channels) {
			plate.mR <- tapply(Rf,plate,hubermu)
			plate.mG <- tapply(Gf,plate,hubermu)
			mR <- mG <- Rf
			for (i in 1:max(plate)) {
				mR[plate==i] <- plate.mR[i]
				mG[plate==i] <- plate.mG[i]
			}
			if(plot) {
				plot(printorder,Rf,xlab="Print Order",ylab="Log Intensity",type="n")
				points(printorder,Rf,pch=".",col="red")
				points(printorder,Gf,pch=".",col="green")
				lines(printorder[ord],mR[ord],col="red")
				lines(printorder[ord],mG[ord],col="green")
				return(invisible())
			}
			mR <- mR - mean(mR,na.rm=TRUE)
			mG <- mG - mean(mG,na.rm=TRUE)
		} else {
			plate.m <- tapply(Af,plate,hubermu)
			m <- Af
			for (i in 1:max(plate)) m[plate==i] <- plate.m[i]
			if(plot) {
				plot(printorder,Af,xlab="Print Order",ylab="Log Intensity",pch=".")
				lines(printorder[ord],m[ord])
			}
			mR <- mG <- m - mean(m,na.rm=TRUE)
		}
	} else {
		# Smooth correction for time order
		if(separate.channels) {
			mR <- fitted(loess(Rf~printorder,span=span,degree=0,family="symmetric",trace.hat="approximate",iterations=5,surface="direct",na.action=na.exclude))
			mG <- fitted(loess(Gf~printorder,span=span,degree=0,family="symmetric",trace.hat="approximate",iterations=5,surface="direct",na.action=na.exclude))
			if(plot) {
				plot(printorder,Rf,xlab="Print Order",ylab="Log Intensity",type="n")
				points(printorder,Rf,pch=".",col="red")
				points(printorder,Gf,pch=".",col="green")
				lines(printorder[ord],mR[ord],col="red")
				lines(printorder[ord],mG[ord],col="green")
				return(invisible())
			}
			mR <- mR - mean(mR,na.rm=TRUE)
			mG <- mG - mean(mG,na.rm=TRUE)
		} else {
			m <- fitted(loess(Af~printorder,span=span,degree=0,family="symmetric",trace.hat="approximate",iterations=5,surface="direct",na.action=na.exclude))
			if(plot) {
				plot(printorder,Af,xlab="Print Order",ylab="Log Intensity",pch=".")
				lines(printorder[ord],m[ord])
			}
			mR <- mG <- m - mean(m,na.rm=TRUE)
		}
	}
	list(R=2^(Rf-mR),G=2^(Gf-mG),R.trend=mR,G.trend=mG)
}

plotPrintorder <- function(object,layout,start="topleft",slide=1,method="loess",separate.channels=FALSE,span=0.1,plate.size=32) {
#	Pre-normalize the foreground intensities for print order.
#	Gordon Smyth
#	9 Apr 2002.  Last revised 18 June 2003.

	if(is.null(object$R) || is.null(object$G)) stop("List must contain components R and G")
	G <- object$G[,slide]
	R <- object$R[,slide]
	if(length(R) != length(G)) stop("R and G must have same length")
	start <- match.arg(start,c("topleft","topright"))
	po <- printorder(layout,start=start)$printorder
	invisible(normalizeForPrintorder.rg(R=R,G=G,printorder=po,method=method,separate.channels=separate.channels,span=span,plate.size=plate.size,plot=TRUE))
}

#  BETWEEN ARRAY NORMALIZATION

normalizeBetweenArrays <- function(object, method="Aquantile", targets=NULL, ...) {
#	Normalize between arrays
#	Gordon Smyth
#	12 Apri 2003.  Last revised 19 October 2005.

	choices <- c("none","scale","quantile","Aquantile","Gquantile","Rquantile","Tquantile","vsn")
	method <- match.arg(method,choices)
	if(method=="vsn") require("vsn")
	if(is(object,"matrix")) {
		if(!(method %in% c("none","scale","quantile","vsn"))) stop("method not applicable to matrix objects")
		return(switch(method,
			none = object,
			scale = normalizeMedianAbsValues(object),
			quantile = normalizeQuantiles(object, ...),
			vsn = vsn(intensities=object,...)@exprs/log(2)
		))
	}
	if(method=="vsn") {
		y <- NULL
		if(!is.null(object$G) && !is.null(object$R)) {
			y <- cbind(object$G,object$R)
			object$G <- object$R <- NULL
		}
		if(!is.null(object$M) && !is.null(object$A)) y <- 2^cbind(object$A-object$M/2,object$A+object$M/2)
		if(is.null(y)) stop("object doesn't appear to be RGList or MAList")
		y <- vsn(intensities=y,...)
		n2 <- ncol(y@exprs)/2
		G <- y@exprs[,1:n2]/log(2)
		R <- y@exprs[,n2+(1:n2)]/log(2)
		object$M <- R-G
		object$A <- (R+G)/2
		if(length(y@description@preprocessing)) object$preprocessing <- y@description@preprocessing
		if(!is(object,"MAList")) object <- new("MAList",unclass(object))
		return(object)
	}
	if(is(object,"RGList")) object <- MA.RG(object)
	if(is.null(object$M) || is.null(object$A)) stop("object doesn't appear to be RGList or MAList object")
	switch(method,
		scale = {
			object$M <- normalizeMedianAbsValues(object$M)
			object$A <- normalizeMedianAbsValues(object$A)
		},
		quantile = {
			narrays <- NCOL(object$M)
			Z <- normalizeQuantiles(cbind(object$A-object$M/2,object$A+object$M/2),...)
			G <- Z[,1:narrays]
			R <- Z[,narrays+(1:narrays)]
			object$M <- R-G
			object$A <- (R+G)/2
		},
		Aquantile = {
			object$A <- normalizeQuantiles(object$A,...)
		},
		Gquantile = {
			G <- object$A-object$M/2
			E <- normalizeQuantiles(G,...) - G
			object$A <- object$A + E
		},
		Rquantile = {
			R <- object$A+object$M/2
			E <- normalizeQuantiles(R,...) - R
			object$A <- object$A + E
		},
		Tquantile = {
			narrays <- NCOL(object$M)
			if(NCOL(targets)>2) targets <- targets[,c("Cy3","Cy5")]
			targets <- as.vector(targets)
			Z <- cbind(object$A-object$M/2,object$A+object$M/2)
			for (u in unique(targets)) {
				j <- targets==u
				Z[,j] <- normalizeQuantiles(Z[,j],...)
			}
			G <- Z[,1:narrays]
			R <- Z[,narrays+(1:narrays)]
			object$M <- R-G
			object$A <- (R+G)/2
		})
	object
}

normalizeQuantiles <- function(A, ties=TRUE) {
#	Normalize columns of a matrix to have the same quantiles, allowing for missing values.
#	Gordon Smyth
#	25 June 2002.  Last revised 8 June 2006.

	n <- dim(A)
	if(is.null(n)) return(A)
	if(n[2]==1) return(A)
	O <- S <- array(,n)
	nobs <- rep(n[1],n[2])
	i <- (0:(n[1]-1))/(n[1]-1)
	for (j in 1:n[2]) {
		Si <- sort(A[,j], method="quick", index.return=TRUE)
		nobsj <- length(Si$x)
		if(nobsj < n[1]) {
			nobs[j] <- nobsj
			isna <- is.na(A[,j])
			S[,j] <- approx((0:(nobsj-1))/(nobsj-1), Si$x, i, ties="ordered")$y
			O[!isna,j] <- ((1:n[1])[!isna])[Si$ix]
		} else {
			S[,j] <- Si$x
			O[,j] <- Si$ix
		}
	}
	m <- rowMeans(S)
	for (j in 1:n[2]) {
		if(ties) r <- rank(A[,j])
		if(nobs[j] < n[1]) {
			isna <- is.na(A[,j])
			if(ties)
				A[!isna,j] <- approx(i, m, (r[!isna]-1)/(nobs[j]-1), ties="ordered")$y
			else
				A[O[!isna,j],j] <- approx(i, m, (0:(nobs[j]-1))/(nobs[j]-1), ties="ordered")$y
		} else {
			if(ties)
				A[,j] <- approx(i, m, (r-1)/(n[1]-1), ties="ordered")$y
			else
				A[O[,j],j] <- m
		}
	}
	A
}

normalizeMedianAbsValues <- function(x) 
{
#	Normalize columns of a matrix to have the same median absolute value
#	Gordon Smyth
#	12 April 2003.  Last modified 19 Oct 2005.

	narrays <- NCOL(x)
	if(narrays==1) return(x)
	cmed <- log(apply(abs(x), 2, median, na.rm=TRUE))
	cmed <- exp(cmed - mean(cmed))
	t(t(x)/cmed)
}


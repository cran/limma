#  SEPARATE CHANNEL ANALYSIS

lmscFit <- function(object,design,correlation)
#	Fit single channel linear model for each gene to a series of microarrays
#	allowing for known correlation between the channels on each spot.
#	Gordon Smyth
#	14 March 2004.  Last modified 26 June 2004.
{
#	Check input
	M <- as.matrix(object$M)
	A <- as.matrix(object$A)
	if(is.null(M) || is.null(A)) stop("object must have components M and A")
	dimM <- dim(M)
	dimA <- dim(A)
	if(any(dimM != dimA)) stop("dimensions of M and A don't match")
	if(!all(is.finite(M)) || !all(is.finite(A))) stop("Missing or infinite values found in M or A")
	if(missing(design)) stop("design matrix must be specified")
	narrays <- dimM[2]
	ny <- 2*narrays
	design <- as.matrix(design)
	if(nrow(design) != ny) stop("The number of rows of the design matrix should match the number of channel intensities, i.e., twice the number of arrays")
	if(missing(correlation)) stop("intra-spot correlation must be specified")
	if(abs(correlation) >= 1) stop("correlation must be strictly between -1 and 1")

#	Dimensions
	nbeta <- ncol(design)
	coef.names <- colnames(design)
	ngenes <- dimM[1]

#	Main computation
	sdM <- sqrt(2*(1-correlation))
	sdA <- sqrt((1+correlation)/2)
	y <- rbind(t(M)/sdM, t(A)/sdA)
	designM <- (diag(narrays) %x% matrix(c(-1,1),1,2)) %*% design
	designA <- (diag(narrays) %x% matrix(c(0.5,0.5),1,2)) %*% design
	X <- rbind(designM/sdM, designA/sdA)

#	In future it may be necessary to allow for quality weights, this call does not
	fit <- lm.fit(X,y)
	fit$sigma <- sqrt(colSums(fit$effects[(fit$rank+1):ny,]^2) / fit$df.residual)
	fit$fitted.values <- fit$residuals <- fit$effects <- NULL
	fit$coefficients <- t(fit$coefficients)
	stdev.unscaled <- sqrt(diag(chol2inv(fit$qr$qr)))
	fit$stdev.unscaled <- matrix(stdev.unscaled,ngenes,nbeta,byrow=TRUE)
	fit$df.residual <- rep.int(fit$df.residual,ngenes)
	dimnames(fit$stdev.unscaled) <- dimnames(fit$stdev.unscaled) <- dimnames(fit$coefficients)
	fit$design <- design
	fit$correlation <- correlation
	fit$genes <- object$genes
	fit$Amean <- rowMeans(A,na.rm=TRUE)
	fit$cov.coefficients <- chol2inv(fit$qr$qr,size=fit$qr$rank)
	fit$pivot <- fit$qr$pivot
	new("MArrayLM",fit)
}

intraspotCorrelation <- function(object,design,trim=0.15)
#	Estimate intra-spot correlation between channels for two channel data
#	Gordon Smyth
#	19 April 2004.  Last modified 13 Nov 2005.
{
#	Check input
	M <- as.matrix(object$M)
	A <- as.matrix(object$A)
	if(is.null(M) || is.null(A)) stop("object should have components M and A")
	dimM <- dim(M)
	dimA <- dim(A)
	if(any(dimM != dimA)) stop("dimensions of M and A don't match")
	if(!all(is.finite(M)) || !all(is.finite(A))) stop("Missing or infinite values found in M or A")
	if(missing(design)) stop("design matrix must be specified")
	ngenes <- dimM[1]
	narrays <- dimM[2]
	ny <- 2*narrays
	design <- as.matrix(design)
	if(nrow(design) != ny) stop("The number of rows of the design matrix should match the number of channel intensities, i.e., twice the number of arrays")

#	Fit heteroscedastic regression for each gene
	Ident <- diag(narrays)
	designM <- (Ident %x% matrix(c(-1,1),1,2)) %*% design
	designA <- (Ident %x% matrix(c(0.5,0.5),1,2)) %*% design
	X <- rbind(designM, designA)
	Z <- diag(2) %x% rep(1,narrays)
	if(!require(statmod)) stop("statmod package required but is not available")
	arho <- rep(NA,ngenes)
	degfre <- matrix(0,ngenes,2,dimnames=list(rownames(M),c("df.M","df.A")))
	for (i in 1:ngenes) {
		y <- c(M[i,],A[i,])
		fit <- try(remlscore(y,X,Z),silent=TRUE)
		if(is.list(fit)) {
			arho[i] <- 0.5*(fit$gamma[2]-fit$gamma[1])
			degfre[i,] <- crossprod(Z,1-fit$h)
		}
	}
	arho <- arho+log(2)
	degfreNA <- degfre
	degfreNA[degfre==0] <- NA
	arhobias <- digamma(degfreNA[,1]/2)-log(degfreNA[,1]/2)-digamma(degfreNA[,2]/2)+log(degfreNA[,2]/2)
	list(consensus.correlation=tanh(mean(arho-arhobias,trim=trim,na.rm=TRUE)), atanh.correlations=arho, df=degfre)
}

targetsA2C <- function(targets,channel.codes=c(1,2),channel.columns=list(Target=c("Cy3","Cy5")),grep=FALSE)
#	Convert data.frame with one row for each two-color array
#	into data.frame with one row for each channel
#	Gordon Smyth
#	16 March 2004.  Last modified 25 May 2004.
{
	targets <- as.data.frame(targets)
	narrays <- nrow(targets)
	nchannelcol <- 0
	nothercol <- ncol(targets)
	if(narrays==0 || nothercol==0) return(targets)
	lcc <- length(channel.columns)
	if(lcc) {
		hyb <- channel.columns
		cheaders <- names(channel.columns)
		for (i in 1:lcc) {
			aheaders <- channel.columns[[i]]
			if(grep) {
				k <- grep(tolower(aheaders[1]),tolower(names(targets)))
				if(length(k)==1) aheaders[1] <- names(targets)[k]
				k <- grep(tolower(aheaders[2]),tolower(names(targets)))
				if(length(k)==1) aheaders[2] <- names(targets)[k]
			}
			if(all(aheaders %in% names(targets))) {
				hyb[[i]] <- as.vector(as.matrix((targets[,aheaders])))
				targets[[ aheaders[1] ]] <- NULL
				targets[[ aheaders[2] ]] <- NULL
				nchannelcol <- nchannelcol+1
				nothercol <- nothercol-2
			} else {
				hyb[[ cheaders[i] ]] <- NULL
			}
		}
	}
	channel.col <- rep(channel.codes,each=narrays)
	out <- data.frame(channel.col=I(channel.col),row.names=paste(row.names(targets),channel.col,sep="."))
	if(nothercol) out <- cbind(out,rbind(targets,targets))
	if(nchannelcol) out <- cbind(out,hyb)
	o <- as.vector(t(matrix(1:(2*narrays),narrays,2)))
	out[o,]
}

designI2M <- function(design)
#  Convert individual channel design matrix to design matrix for log-ratios
#  Gordon Smyth
#  22 June 2004
{
	design <- as.matrix(design)
	narrays <- nrow(design)/2
	(diag(narrays) %x% matrix(c(-1,1),1,2)) %*% design
}

designI2A <- function(design)
#  Convert individual channel design matrix to design matrix for A-values
#  Gordon Smyth
#  22 June 2004
{
	design <- as.matrix(design)
	narrays <- nrow(design)/2
	(diag(narrays) %x% matrix(c(0.5,0.5),1,2)) %*% design
}

exprs.MA <- function(MA)
#	Extract matrix of log-expression data from MAList object
#	Gordon Smyth, 29 August 2006
{
	y <- rbind(as.matrix(MA$A-MA$M/2),as.matrix(MA$A+MA$M/2))
	dim(y) <- c(nrow(y)/2,ncol(y)*2)
	y
}

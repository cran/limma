#  EVALUATION

bwss <- function(x,group) {
#	Within and between group sums of squares
#	Gordon Smyth
#	14 Mar 2002. Revised 13 May 2002.

	is.na.xg <- is.na(x) | is.na(group)
	if(any(is.na.xg)) {
		x <- x[!is.na.xg]
		group <- group[!is.na.xg]
	}
	if(length(x)==0) return(list(bss=NA,wss=NA,bdf=NA,wdf=NA))
	m <- tapply(x,group,mean)
	v <- tapply(x,group,var)
	n <- table(group)
	if(any(n==0)) {
		m <- m[n>0]
		v <- v[n>0]
		n <- n[n>0]
	}
	mm <- sum(n*m)/sum(n)
	wss <- sum((n-1)*v,na.rm=TRUE)
	bss <- sum(n*(m-mm)^2)
	wdf <- sum(n-1)
	bdf <- length(m)-1
	list(bss=bss,wss=wss,bdf=bdf,wdf=wdf)
}

bwss.matrix <- function(x) {
#	Sums of squares between and within columns of a matrix
#	Gordon Smyth
#	16 Mar 2002. Revised 13 May 2002.

	n <- colSums(!is.na(x))
	if(any(n==0)) {
		x <- x[,n>0,drop=FALSE]
		n <- n[n>0]
	}
	nc <- length(n)
	if(nc==0) return(list(bss=NA,wss=NA,bdf=NA,wdf=NA))
	m <- colMeans(x,na.rm=TRUE)
	v <- apply(x,2,var,na.rm=TRUE)
	mm <- sum(n*m)/sum(n)
	wss <- sum((n-1)*v,na.rm=TRUE)
	bss <- sum(n*(m-mm)^2)
	wdf <- sum(n-1)
	bdf <- nc-1
	list(bss=bss,wss=wss,bdf=bdf,wdf=wdf)
}

anova.MAList <- function(object,design=NULL,ndups=2,...) {
#	Analysis of variance, gene x array, for series of replicate arrays
#	Useful for comparing between to within gene variability
#	Gordon Smyth
#	16 March 2002.  Last modified 17 Aug 2004.

	M <- object$M
	if(!is.null(design)) {
		d <- dim(design)
		if(!is.null(d)) if(d[2]>1) stop("Design matrix should have only one column")
		M <- .matvec(M,design)
	}
	bwss.array <- bwss.matrix(M)
	nspots <- dim(M)[1]
	narrays <- dim(M)[2]
	ngenes <- nspots/ndups
	dim(M) <- c(ndups,ngenes*narrays)
	bwss.genearray <- bwss.matrix(M)
	dim(M) <- c(ndups,ngenes,narrays)
	M <- aperm(M,c(1,3,2))
	dim(M) <- c(ndups*narrays,ngenes)
	bwss.gene <- bwss.matrix(M)
	src <- c("Genes","Arrays","Interaction","Duplicates")
	ss <- c(bwss.gene$bss,bwss.array$bss,bwss.genearray$bss-bwss.array$bss-bwss.gene$bss,bwss.genearray$wss)
	df <- c(bwss.gene$bdf,bwss.array$bdf,bwss.genearray$bdf-bwss.array$bdf-bwss.gene$bdf,bwss.genearray$wdf)
	ms <- ss/df
	sig <- ms
	sig[1] <- (ms[1]-ms[3])/ndups/narrays
	sig[2] <- (ms[2]-ms[3])/ndups/ngenes
	sig[3] <- (ms[3]-ms[4])/ndups
	ratio <- c(sig[1]/sum(sig[2:4]),NA,NA,NA)
	table <- data.frame(df, ss, ms, sig, ratio)
	dimnames(table) <- list(src,c("Df","Sum Sq","Mean Sq","Var Comp","Ratio")) 
	structure(table,heading=c("Analysis of Variance Table\n","Between and Within Genes"),class=c("anova","data.frame"))
}

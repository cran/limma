#  TOPTABLE.R

topTable <- function(fit,coef=NULL,number=10,genelist=fit$genes,adjust.method="BH",sort.by="B",resort.by=NULL)
#	Summary table of top genes, object-orientated version
#	Gordon Smyth
#	4 August 2003.  Last modified 27 Oct 2006.
{
	if(is.null(coef)) {
		if(ncol(fit)>1)
			return(topTableF(fit,number=number,genelist=genelist,adjust.method=adjust.method))
		else
			coef=1
	}
	if(length(coef)>1) return(topTableF(eBayes(fit[,coef]),number=number,genelist=genelist,adjust.method=adjust.method))
	fit <- unclass(fit)
	toptable(fit=fit[c("coefficients","stdev.unscaled")],
		coef=coef,
		number=number,
		genelist=genelist,
		A=fit$Amean,
		eb=fit[c("t","p.value","lods")],
		adjust.method=adjust.method,
		sort.by=sort.by,
		resort.by=resort.by)
}

topTableF <- function(fit,number=10,genelist=fit$genes,adjust.method="BH")
#	Summary table of top genes by F-statistic
#	Gordon Smyth
#	27 August 2006. Last modified 27 Oct 2006.
{
#	Check input
	if(!is.null(genelist) && is.null(dim(genelist))) genelist <- data.frame(ProbeID=I(genelist))
	M <- as.matrix(fit$coefficients)
	if(is.null(colnames(M))) colnames(M) <- paste("Coef",1:ncol(M),sep="")

	adj.P.Value <- p.adjust(fit$F.p.value,method=adjust.method)
	o <- order(fit$F.p.value,decreasing=FALSE)[1:number]
	if(is.null(genelist))
		tab <- data.frame(M[o,,drop=FALSE])
	else
		tab <- data.frame(genelist[o,,drop=FALSE],M[o,,drop=FALSE])
	if(!is.null(fit$Amean)) tab <- data.frame(tab,AveExpr=fit$Amean[o])
	tab <- data.frame(tab,F=fit$F[o],P.Value=fit$F.p.value[o],adj.P.Val=adj.P.Value[o])
	rownames(tab) <- as.character(1:length(M))[o]
	tab
}

toptable <- function(fit,coef=1,number=10,genelist=NULL,A=NULL,eb=NULL,adjust.method="BH",sort.by="B",resort.by=NULL,...)
#	Summary table of top genes
#	Gordon Smyth
#	21 Nov 2002. Last revised 27 Mar 2007.
{
	if(is.null(eb)) {
		fit$coefficients <- as.matrix(fit$coefficients)[,coef]
		fit$stdev.unscaled <- as.matrix(fit$stdev.unscaled)[,coef]
		eb <- ebayes(fit,...)
		coef <- 1
	}
	M <- as.matrix(fit$coefficients)[,coef]
	if(length(M) < number) number <- length(M)
	if(number < 1) return(data.frame())
	if(is.null(A)) {
		if(sort.by=="A") stop("Cannot sort by A-values as these have not been given")
	} else {
		if(NCOL(A)>1) A <- rowMeans(A,na.rm=TRUE)
	}
	tstat <- as.matrix(eb$t)[,coef]
	P.Value <- as.matrix(eb$p.value)[,coef]
	B <- as.matrix(eb$lods)[,coef]
	sort.by <- match.arg(sort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t","B"))
	if(sort.by=="M") sort.by="logFC"
	if(sort.by=="A" || sort.by=="Amean") sort.by <- "AveExpr"
	ord <- switch(sort.by,
		logFC=order(abs(M),decreasing=TRUE),
		AveExpr=order(A,decreasing=TRUE),
		P=order(P.Value,decreasing=FALSE),
		p=order(P.Value,decreasing=FALSE),
		T=order(abs(tstat),decreasing=TRUE),
		t=order(abs(tstat),decreasing=TRUE),
		B=order(B,decreasing=TRUE)
	)
	top <- ord[1:number]
	adj.P.Value <- p.adjust(P.Value,method=adjust.method)
	if(is.null(genelist))
		tab <- data.frame(logFC=M[top])
	else {
		if(is.null(dim(genelist)))
			tab <- data.frame(ID=I(genelist[top]),logFC=M[top])
		else
			tab <- data.frame(genelist[top,,drop=FALSE],logFC=M[top])
	}
	if(!is.null(A)) tab <- data.frame(tab,AveExpr=A[top])
	tab <- data.frame(tab,t=tstat[top],P.Value=P.Value[top],adj.P.Val=adj.P.Value[top],B=B[top])
	rownames(tab) <- as.character(1:length(M))[top]
	if(!is.null(resort.by)) {
		resort.by <- match.arg(resort.by,c("logFC","M","A","Amean","AveExpr","P","p","T","t","B"))
		if(resort.by=="M") resort.by <- "logFC"
		if(resort.by=="A" || resort.by=="Amean") resort.by <- "AveExpr"
		if(resort.by=="p") resort.by <- "P"
		if(resort.by=="T") resort.by <- "t"
		ord <- switch(resort.by,
			logFC=order(tab$logFC,decreasing=TRUE),
			AveExpr=order(tab$AveExpr,decreasing=TRUE),
			P=order(tab$P.Value,decreasing=FALSE),
			t=order(tab$t,decreasing=TRUE),
			B=order(tab$B,decreasing=TRUE)
		)
		tab <- tab[ord,]
	}
#	attr(tab,"coef") <- coef
#	attr(tab,"adjust.method") <- adjust.method
	tab
}


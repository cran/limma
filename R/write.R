#  OUTPUT

write.fit <- function(fit, results=NULL, FStat=NULL, file, digits=2, adjust="none") {
#	Write an MArrayLM fit to a file
#	Gordon Smyth
#	14 Nov 2003.  Last modified 6 July 2004.

	if(!is(fit, "MArrayLM")) stop("fit should be a MArrayLM object")
	if(!is.null(results) && !is(results,"TestResults")) stop("results should be a TestResults object")
	if(is.null(fit$t) || is.null(fit$p.value)) fit <- eBayes(fit)
	p.value <- as.matrix(fit$p.value)
	for (j in 1:ncol(p.value)) p.value[,j] <- p.adjust(p.value[,j],method=adjust)
	tab <- data.frame(A=round(fit$Amean,digits=digits-1),round(fit$coef,digits=digits),round(fit$t,digits=digits-1),round(p.value,digits=digits+2),check.names=FALSE)
	if(!is.null(FStat)) tab <- data.frame(tab,F=round(FStat,digits=digits-1),check.names=FALSE)
	if(!is.null(results)) tab <- data.frame(tab,unclass(results),check.names=FALSE)
	tab <- data.frame(tab,fit$genes,check.names=FALSE)
	write.table(tab,file=file,quote=FALSE,row.names=FALSE,sep="\t")
}

#  OUTPUT

write.fit <- function(fit, results=NULL, file, digits=2, adjust="none", sep="\t") {
#	Write an MArrayLM fit to a file
#	Gordon Smyth
#	14 Nov 2003.  Last modified 24 August 2005.

	if(!is(fit, "MArrayLM")) stop("fit should be an MArrayLM object")
	if(!is.null(results) && !is(results,"TestResults")) stop("results should be a TestResults object")

	if(is.null(fit$t) || is.null(fit$p.value)) fit <- eBayes(fit)
	p.value <- as.matrix(fit$p.value)
	for (j in 1:ncol(p.value)) p.value[,j] <- p.adjust(p.value[,j],method=adjust)

	rn <- function(x,digits=digits) if(is.null(x)) NULL else round(x,digits=digits)
	tab <- list()
	tab$A <- rn(fit$Amean,digits=digits-1)
	tab$Coef <- rn(fit$coef,digits=digits)
	tab$t <- rn(fit$t,digits=digits-1)
	tab$p.value <- rn(p.value,digits=digits+2)
	tab$F <- rn(fit$F,digits=digits-1)
	tab$F.p.value <- rn(fit$F.p.value,digits=digits+2)
	tab$Res <- unclass(results)
	tab$Genes <- fit$genes
	tab <- data.frame(tab,check.names=FALSE)
	write.table(tab,file=file,quote=FALSE,row.names=FALSE,sep=sep)
}

as.data.frame.MArrayLM <- function(x, row.names = NULL, optional = FALSE)
#	Convert MAList object to data.frame
#	Gordon Smyth
#	6 April 2005.  Last modified 24 August 2005.
{
	x <- unclass(x)
	if(is.null(x$coefficients)) {
		warning("NULL coefficients, returning empty data.frame")
		return(data.frame())
	}
	cn <- names(x)
	nprobes <- NROW(x$coefficients)
	include.comp <- cn[unlist(lapply(x,NROW))==nprobes]
	other.comp <- setdiff(names(x),include.comp)
	if(length(other.comp)) for (a in other.comp) x[[a]] <- NULL
#	coef.comp <- c("coefficients","stdev.unscaled","t","p.value","lods")
#	for (a in coef.comp) if(!is.null(x[[a]]) && NCOL(x[[a]])==1) colnames(x[[a]]) <- paste(a,colnames(x[[a]]),sep=".")
	x <- lapply(x,drop)
	as.data.frame(x,row.names=row.names,optional=optional)
}

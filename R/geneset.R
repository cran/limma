##  GENESET.R

geneSetTest <- function(selected,statistics,alternative="mixed",type="auto",ranks.only=TRUE,nsim=10000)
#	Gene set test using either Wilcox test or simulation.
#	Gordon Smyth
#	3 September 2004. Last modified 21 July 2006.
{
	alternative <- match.arg(alternative,c("mixed","either","down","up","less","greater"))
	if(alternative=="two.sided") alternative <- "either"
	if(alternative=="less") alternative <- "down"
	if(alternative=="greater") alternative <- "up"
	type <- match.arg(tolower(type),c("auto","t","f"))
	allsamesign <- all(statistics >= 0) || all(statistics <= 0)
	if(type=="auto") {
		if(allsamesign)
			type <- "f"
		else
			type <- "t"
	}
	if(type=="f" & alternative!="mixed") stop("Only alternative=\"mixed\" is possible with F-like statistics.")
	if(alternative=="mixed") statistics <- abs(statistics)
	if(alternative=="down") {
		statistics <- -statistics
		alternative <- "up"
	}
	if(ranks.only) {
#		The test statistic is the mean rank of the selected statistics
#		and the p-value is obtained explicitly from the Wilcox test
		if(alternative=="either")
			wilc.alt <- "two.sided"
		else
			wilc.alt <- "greater"
		x <- y <- NULL
		if(is.logical(selected)) {
			x <- statistics[selected]
			y <- statistics[!selected]
		}
		if(is.numeric(selected)) {
			x <- statistics[selected]
			y <- statistics[-selected]
		}
		if(is.character(selected)) {
			nam <- names(statistics)
			if(is.null(nam)) stop("selected is character but elements of statistics are not named")
			selected <- is.element(nam,selected)
			x <- statistics[selected]
			y <- statistics[!selected]
		}
		return(wilcox.test(x,y,alternative=wilc.alt,conf.int=FALSE)$p.value)
	} else {
#		The test statistic is the mean of the selected statistics
#		and the p-value is obtained by random permutation
		ssel <- statistics[selected]
		ssel <- ssel[!is.na(ssel)]
		nsel <- length(ssel)
		if(nsel==0) return(1)
		stat <- statistics[!is.na(statistics)]
		msel <- mean(ssel)
		if(alternative=="either")
			posstat <- abs
		else
			posstat <- function(x) x
		msel <- posstat(msel)
		ntail <- 0
		for (i in 1:nsim) if(posstat(mean(sample(stat,nsel))) >= msel) ntail <- ntail+1
		return(ntail/nsim)
	}
}

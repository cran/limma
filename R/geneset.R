##  GENESET.R

geneSetTest <- function(selected,statistics,alternative="two.sided",nsim=10000,ranks.only=FALSE)
#	Gene set test
#	Gordon Smyth
#	3 September 2004
#	Last modified 18 Sep 2004.
{
	alternative <- match.arg(alternative,c("two.sided","less","greater"))
	if(ranks.only) {
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
		return(wilcox.test(x,y,alternative=alternative,conf.int=FALSE)$p.value)
	} else {
		ssel <- statistics[selected]
		ssel <- ssel[!is.na(ssel)]
		nsel <- length(ssel)
		if(nsel==0) return(1)
		stat <- statistics[!is.na(statistics)]
		msel <- mean(statistics[selected],na.rm=TRUE)
		posstat <- switch(alternative,
			"two.sided"=abs,
			"less"=function(x) -x,
			"greater"=function(x) x
		)
		msel <- posstat(msel)
		ntail <- 0
		for (i in 1:nsim) if(posstat(mean(sample(stat,nsel))) >= msel) ntail <- ntail+1
		return(ntail/nsim)
	}
}

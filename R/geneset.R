##  GENESET.R

geneSetTest <- function(selected,statistics,nsim=10000)
#	Gene set test
#	Gordon Smyth
#	3 September 2004
{
	ssel <- statistics[selected]
	ssel <- ssel[!is.na(ssel)]
	nsel <- length(ssel)
	if(nsel==0) return(1)
	stat <- statistics[!is.na(statistics)]
	msel <- mean(statistics[selected])
	ntail <- 0
	for (i in 1:nsim) ntail <- ntail+(mean(sample(stat,nsel))>=msel)
	ntail/nsim
}

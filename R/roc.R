#  ROC.R

auROC <- function(truth, stat) {
#	Area under Receiver Operating Curve for empirical data
#	Gordon Smyth
#	21 Dec 2003

	if(!length(truth)) return(NULL)
    if(!all(sort(unique(truth)) == c(0, 1)))
        stop("'truth' must take values 0 or 1")
	if(!missing(stat)) {
		if(length(stat) != length(truth)) stop("lengths differ")
		isna <- is.na(stat)
		if(any(isna)) truth[isna] <- NA 
		truth <- truth[order(stat,decreasing=TRUE)]
	}
	isna <- is.na(truth)
	if(any(isna)) truth <- truth[!isna] 
	sens <- cumsum(truth)/sum(truth)
	mean(sens[truth==0])
}


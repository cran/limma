#  POOL VAR

poolVar <- function(var, df=n-1, multiplier=1/n, n)
#	Pool sample variances for unequal variances
#	as for Welch's test
#	Gordon Smyth
#	22 Jan 2004.
{
	s2 <- as.vector(var)
	df <- as.vector(df)
	m <- as.vector(multiplier)
	if(any(m<0)) stop("Multipliers must be non-negative")
	if(all(m==0)) return(list(var=0,df=0,multiplier=0))
	sm <- sum(m)
	m <- m/sm
	out <- list(var = sum(m * s2))
	out$df <- out$var^2 / sum(m^2 * s2^2 / df)
	out$multiplier <- sm
	out
}


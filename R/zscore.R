#  SCORE.R

zscoreGamma <- function(q, shape, rate = 1, scale = 1/rate) 
#  Z-score equivalents for gamma deviates
#  Gordon Smyth
#  1 October 2003
{
	z <- q
	n <- length(q)
	shape <- rep(shape,length.out=n)
	scale <- rep(scale,length.out=n)
	up <- (q > shape*scale)
	if(any(up)) z[up] <- qnorm(pgamma(q[up],shape=shape[up],scale=scale[up],lower.tail=FALSE,log.p=TRUE),,lower.tail=FALSE,log.p=TRUE)
	if(any(!up)) z[!up] <- qnorm(pgamma(q[!up],shape=shape[!up],scale=scale[!up],lower.tail=TRUE,log.p=TRUE),lower.tail=TRUE,log.p=TRUE)
	z
}

zscoreT <- function(x, df)
#  Z-score equivalents for gamma deviates
#  Gordon Smyth
#  24 August 2003
{
	z <- x
	df <- rep(df,length.out=length(x))
	pos <- x>0
	if(any(pos)) z[pos] <- qnorm(pt(x[pos],df=df[pos],lower.tail=FALSE,log.p=TRUE),lower.tail=FALSE,log.p=TRUE) 
	if(any(!pos)) z[!pos] <- qnorm(pt(x[!pos],df=df[!pos],lower.tail=TRUE,log.p=TRUE),lower.tail=TRUE,log.p=TRUE)
	z
}

tZscore <- function(x, df)
#  t-statistics equivalents for z-scores deviates
#  Gordon Smyth
#  1 June 2004
{
	z <- x
	df <- rep(df,length.out=length(x))
	pos <- x>0
	if(any(pos)) z[pos] <- qt(pnorm(x[pos],lower.tail=FALSE,log.p=TRUE),df=df[pos],lower.tail=FALSE,log.p=TRUE) 
	if(any(!pos)) z[!pos] <- qt(pnorm(x[!pos],lower.tail=TRUE,log.p=TRUE),df=df[!pos],lower.tail=TRUE,log.p=TRUE)
	z
}

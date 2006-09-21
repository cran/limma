#  plotlines.R

plotlines <- function(x,first.column.origin=FALSE,xlab="Column",ylab="x",col="black",lwd=1,...)
#	Time-course-style plot
#	Plot expression data as a line for each probe
#	Gordon Smyth
#	27 August 2006. Last revised 30 August 2006.
{
	x <- as.matrix(x)
	ntime <- ncol(x)
	ngenes <- nrow(x)
	if(first.column.origin) x <- x - array(x[,1],dim(x))
	time <- matrix(1:ntime,ngenes,ntime,byrow=TRUE)
	plot(time,x,type="n",xlab=xlab,ylab=ylab,...)
	col <- rep(col,ngenes)
	for (i in 1:ngenes) lines(1:ntime,x[i,],col=col[i],lwd=lwd)
}

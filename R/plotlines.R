plotlines <- function(x,first.column.origin=FALSE,xlab="Column",ylab="x",...)
#	Time-course-style plot
#	Plot expression data as a line for each probe
#	Gordon Smyth, 27 August 2006
{
	x <- as.matrix(x)
	ntime <- ncol(x)
	ngenes <- nrow(x)
	if(normalize) x <- x - array(x[,1],dim(x))
	time <- matrix(1:ntime,ngenes,ntime,byrow=TRUE)
	plot(time,x,type="n",xlab=xlab,ylab=ylab,...)
	for (i in 1:ngenes) lines(1:ntime,x[i,])
}

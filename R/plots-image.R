#  PRESENTATION PLOTS

imageplot <- function(z, layout, low=NULL, high=NULL, ncolors=123, zerocenter=NULL, zlim=NULL, mar=c(2,1,1,1), legend=TRUE, ...) {
#  Image plot of spotted microarray data
#  Gordon Smyth
#  20 Nov 2001.  Last revised 18 Sep 2004.

#  Check input
	gr <- layout$ngrid.r
	gc <- layout$ngrid.c
	sr <- layout$nspot.r
	sc <- layout$nspot.c
	if(is.null(gr)||is.null(gc)||is.null(sr)||is.null(sc)) stop("Layout needs to contain components ngrid.r, ngrid.c, nspot.r and nspot.c")
	if(length(z) != gr*gc*sr*sc) stop("Number of image spots does not agree with layout dimensions")

#  Check colours
	if(is.character(low)) low <- col2rgb(low)/255
	if(is.character(high)) high <- col2rgb(high)/255
	if(!is.null(low) && is.null(high)) high <- c(1,1,1) - low
	if(is.null(low) && !is.null(high)) low <- c(1,1,1) - high

#  Is zlim preset?
	zr0 <- range(z,na.rm=TRUE)
	if(!is.null(zlim)) {
		z <- pmax(zlim[1],z)
		z <- pmin(zlim[2],z)
	}

#  Plot differential expression from "green" to "red" or plot one variable from "white" to "blue"?
	zr <- range(z,na.rm=TRUE)
	if(!all(is.finite(zr))) stop("Infinite values found: consider using finite zlim value")
	zmax <- max(abs(zr))
	zmin <- zr[1]
	if(is.null(zerocenter)) zerocenter <- (zmin < 0)
	if(zerocenter) {
		if(is.null(low)) low <- c(0,1,0)
		if(is.null(high)) high <- c(1,0,0)
		if(is.null(zlim)) zlim <- c(-zmax,zmax)
	} else {
		if(is.null(low)) low <- c(1,1,1)
		if(is.null(high)) high <- c(0,0,1)
		if(is.null(zlim)) zlim <- c(zmin,zmax)
	}

#  Now make the plot
	col <- rgb( seq(low[1],high[1],len=ncolors), seq(low[2],high[2],len=ncolors), seq(low[3],high[3],len=ncolors) )
	dim(z) <- c(sc,sr,gc,gr)
	z <- aperm(z,perm=c(2,4,1,3))
	dim(z) <- c(gr*sr,gc*sc)
	old.par <- par(mar=mar)
	on.exit(par(old.par))
	image(0:(gr*sr),0:(gc*sc),z,zlim=zlim,col=col,xaxt="n",yaxt="n",...)
	for (igrid in 0:gc) lines( c(0,gr*sr), rep(igrid*sc,2) )
	for (igrid in 0:gr) lines( rep(igrid*sr,2), c(0,gc*sc) )
	if(legend) mtext(paste("z-range ",round(zr0[1],1)," to ",round(zr0[2],1)," (saturation ",round(zlim[1],1),", ",round(zlim[2],1),")",sep=""),side=1,cex=0.6)
	invisible()
}

imageplot3by2 <- function(RG, z="Gb", prefix=paste("image",z,sep="-"), path=NULL, zlim=NULL, common.lim=TRUE, ...)
#	Make files of image plots, six to a page
#	Gordon Smyth  10 June 2004.
#	Suggestions of Marcus Davy implemented 30 September 2005.
{
	if(is.null(path)) path="."
	narrays <- ncol(RG)
	npages <- ceiling(narrays/6)
	cnames <- colnames(RG)
	if(is.null(zlim) && common.lim) zlim <- quantile(RG[[z]],c(0.05,0.95),na.rm=TRUE)
	for (ipage in 1:npages) {
		i1 <- ipage*6-5
		i2 <- min(ipage*6,narrays)
		fn <- file.path(path, paste(prefix,"-",i1,"-",i2,".png",sep=""))
		png(filename=fn,width=6.5*140,height=10*140,pointsize=20)
		par(mfrow=c(3,2))
		for (i in i1:i2) {
			imageplot(RG[[z]][,i],RG$printer,zlim=zlim,mar=c(2,2,4,2),main=cnames[i],...)
		}
		dev.off()
	}
	invisible()
}


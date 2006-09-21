#  VENN DIAGRAM COUNTS AND PLOTS

vennCounts <- function(x,include="both") {
#	Venn diagram counts
#	Gordon Smyth  
#	4 July 2003.  Last modified 30 May 2004.

	x <- as.matrix(x)
	include <- match.arg(include,c("both","up","down"))
	x <- sign(switch(include,
		both = abs(x),
		up = x > 0,
		down = x < 0
	))
	nprobes <- nrow(x)
	ncontrasts <- ncol(x)
	names <- colnames(x)
	if(is.null(names)) names <- paste("Group",1:ncontrasts)
	noutcomes <- 2^ncontrasts
	outcomes <- matrix(0,noutcomes,ncontrasts)
	colnames(outcomes) <- names
	for (j in 1:ncontrasts)
		outcomes[,j] <- rep(0:1,times=2^(j-1),each=2^(ncontrasts-j))
	xlist <- list()
	for (i in 1:ncontrasts) xlist[[i]] <- factor(x[,ncontrasts-i+1],levels=c(0,1))
	counts <- as.vector(table(xlist))
	structure(cbind(outcomes,Counts=counts),class="VennCounts")
}

vennDiagram <- function(object,include="both",names,mar=rep(1,4),cex=1.5,lwd=1,circle.col,counts.col,show.include,...)
#	Plot Venn diagram
#	Gordon Smyth, James Wettenhall.
#	Capabilities for multiple counts and colors by Francois Pepin.
#	4 July 2003.  Last modified 21 September 2006.
{
	if (!is(object, "VennCounts")){
		if (length(include)>2) stop("Cannot plot Venn diagram for more than 2 sets of counts")
		if (length(include)==2) object.2 <- vennCounts(object, include = include[2])
		object <- vennCounts(object, include = include[1])
	}
	else if(length(include==2)) include <- include[1]
	nsets <- ncol(object)-1
	if(nsets > 3) stop("Can't plot Venn diagram for more than 3 sets")
	if(missing(names)) names <- colnames(object)[1:nsets]
	counts <- object[,"Counts"]
	if(length(include)==2) counts.2 <- object.2[, "Counts"]
	if(missing(circle.col)) circle.col <- par('col')
	if(length(circle.col)<nsets) circle.col <- rep(circle.col,length.out=nsets)
	if(missing(counts.col)) counts.col <- par('col')
	if(length(counts.col)<length(include)) counts.col <- rep(counts.col,length.out=length(include))
	if(missing(show.include)) show.include <- as.logical(length(include)-1)
	theta <- 2*pi*(1:360)/360
	xcentres <- list(0,c(-1,1),c(-1,1,0))[[nsets]]
	ycentres <- list(0,c(0,0),c(1/sqrt(3),1/sqrt(3),-2/sqrt(3)))[[nsets]]
	r <- c(1.5,1.5,1.5)[nsets]
	xtext <- list(-1.2,c(-1.2,1.2),c(-1.2,1.2,0))[[nsets]]
	ytext <- list(1.8,c(1.8,1.8),c(2.4,2.4,-3))[[nsets]]
	old.par <- par(mar=mar)
	on.exit(par(old.par))
	plot(x=0,y=0,type="n",xlim=c(-4,4),ylim=c(-4,4),xlab="",ylab="",axes=FALSE,...)
	for(circle in 1:nsets) {
#		lines() better than symbols() to make circles follow aspect ratio of plot
		lines(xcentres[circle]+r*cos(theta),ycentres[circle]+r*sin(theta),lwd=lwd,col=circle.col[circle])
		text(xtext[circle],ytext[circle],names[circle],cex=cex)
	}
	switch(nsets,
		{
			rect(-3,-2.5,3,2.5)
			printing <- function(counts, cex, adj,col,leg){
				text(2.3,-2.1,counts[1],cex=cex,col=col,adj=adj)
				text(0,0,counts[2],cex=cex,col=col,adj=adj)
				if(show.include) text(-2.3,-2.1,leg,cex=cex,col=col,adj=adj)
			}
			
		}, {
			rect(-3,-2.5,3,2.5)
			printing <- function(counts, cex, adj,col,leg){
				text(2.3,-2.1,counts[1],cex=cex,col=col,adj=adj)
				text(1.5,0.1,counts[2],cex=cex,col=col,adj=adj)
				text(-1.5,0.1,counts[3],cex=cex,col=col,adj=adj)
				text(0,0.1,counts[4],cex=cex,col=col,adj=adj)
				if(show.include) text(-2.3,-2.1,leg,cex=cex,col=col,adj=adj)
			}
		}, {
			rect(-3,-3.5,3,3.3)
			printing <- function(counts, cex, adj,col,leg){
				text(2.5,-3,counts[1],cex=cex,col=col,adj=adj)
				text(0,-1.7,counts[2],cex=cex,col=col,adj=adj)
				text(1.5,1,counts[3],cex=cex,col=col,adj=adj)
				text(.75,-.35,counts[4],cex=cex,col=col,adj=adj)
				text(-1.5,1,counts[5],cex=cex,col=col,adj=adj)
				text(-.75,-.35,counts[6],cex=cex,col=col,adj=adj)
				text(0,.9,counts[7],cex=cex,col=col,adj=adj)
				text(0,0,counts[8],cex=cex,col=col,adj=adj)
				if(show.include) text(-2.5,-3,leg,cex=cex,col=col,adj=adj)
			}
		}
	)
	adj <- c(0.5,0.5)
	if (length(include)==2)
	  adj <- c(0.5,0)
	printing(counts,cex,adj,counts.col[1],include[1])
	if (length(include)==2) printing(counts.2,cex,c(0.5,1),counts.col[2],include[2])
	invisible()
}


#  PRESENTATION PLOTS

heatDiagram <- function(results,coef,primary=1,names=NULL,treatments=colnames(coef),limit=NULL,orientation="landscape",cex=1,low="green",high="red",ncolors=123,...) {
#	Heat diagram to display fold changes of genes under different conditions
#	Gordon Smyth
#	27 Oct 2002. Last revised 18 Sep 2004.

#	Check input
	results <- as.matrix(results)
	results[is.na(results)] <- 0
	coef <- as.matrix(coef)
	if(!identical(dim(results),dim(coef))) stop("results and coef must be the same size")
	nt <- ncol(results)
	if(is.null(names)) names <- as.character(1:nrow(coef))
	names <- substring(names,1,15)
	if(is.null(treatments)) treatments <- as.character(1:nt)
	orientation <- match.arg(orientation,c("landscape","portrait"))

#	Sort coefficients
	DE <- (abs(results[,primary]) > 0.5)
	ng <- sum(DE)
	if(ng == 0) {
		warning("Nothing significant to plot")
		return(invisible())
	}
	results <- results[DE,,drop=FALSE]
	coef <- coef[DE,,drop=FALSE]
	coef[abs(results) < 0.5] <- NA
	names <- names[DE]
	ord <- order(coef[,primary],decreasing=TRUE)

#	Truncate coefficients if limit is preset
	if(!is.null(limit))
		if(limit > 0) {
			coef[coef < -limit] <- -limit
			coef[coef > limit] <- limit
		} else
			warning("limit ignored because not positive")

#	Check colours
	if(is.character(low)) low <- col2rgb(low)/255
	if(is.character(high)) high <- col2rgb(high)/255
	col <- rgb( seq(low[1],high[1],len=ncolors), seq(low[2],high[2],len=ncolors), seq(low[3],high[3],len=ncolors) )

#	Heat plot
	coef <- coef[ord,,drop=FALSE]
	names <- names[ord]
	out <- data.frame(Name=names,coef)
#	Insert white space between up and down
	nup <- sum(coef[,primary]>=0)
	coef <- rbind(coef[1:nup,,drop=FALSE],matrix(NA,1,nt),coef[(nup+1):ng,,drop=FALSE])
	names <- c(names[1:nup],"",names[(nup+1):ng])
	if(orientation=="portrait") {
		coef <- t(coef)
		coef <- coef[,ng:1,drop=FALSE]
	}
	on.exit(par(old.par))
	if(orientation=="portrait") {
		old.par <- par(mar=cex*c(1,1,4,3))
		image(coef,col=col,xaxt="n",yaxt="n",...)
		cext <- cex*min(1,8/nt)
		mtext(paste(" ",treatments,sep=""),side=3,las=2,at=(cext-1)*0.005+(0:(nt-1))/(nt-1),cex=cext)
		cex <- cex*min(1,40/ng)
		mtext(paste(" ",names,sep=""),side=4,las=2,at=(1-cex)*0.005+((ng-1):0)/(ng-1),cex=cex)
	} else {
		old.par <- par(mar=cex*c(5,6,1,1))
		image(coef,col=col,xaxt="n",yaxt="n",...)
		cext <- cex*min(1,12/nt)
		mtext(paste(treatments," ",sep=""),side=2,las=1,at=(1-cext)*0.005+(0:(nt-1))/(nt-1),cex=cext)
		cex <- cex*min(1,60/ng)
		mtext(paste(names," ",sep=""),side=1,las=2,at=(cex-1)*0.005+(0:(ng-1))/(ng-1),cex=cex)
	}
	invisible(out)
}

heatdiagram <- function(stat,coef,primary=1,names=NULL,treatments=colnames(stat),critical.primary=4,critical.other=3,limit=NULL,orientation="landscape",cex=1,low="green",high="red",ncolors=123,...) {
#	Heat diagram to display fold changes of genes under different conditions
#	Similar to heatDiagram with classifyTestsT(stat,t1=critical.primary,t2=critical.other)
#	except that heatdiagram() requires primary column to be significant at the first step-down level
#	Gordon Smyth
#	27 Oct 2002. Last revised 25 Feb 2003.

#	Check input
	stat <- as.matrix(stat)
	coef <- as.matrix(coef)
	if(any(dim(stat) != dim(coef))) stop("stat and coef must be the same size")
	nt <- ncol(stat)
	if(is.null(names)) names <- as.character(1:nrow(stat))
	names <- substring(names,1,15)
	if(is.null(treatments)) treatments <- as.character(1:nt)

#  Sort coefficients
	DE <- (stat[,primary] > critical.primary)
	if(any(is.na(DE))) DE[is.na(DE)] <- FALSE
	ng <- sum(DE)
	if(sum(DE) == 0) {
		warning("Nothing significant to plot")
		return(invisible())
	}
	stat <- stat[DE,,drop=FALSE]
	coef <- coef[DE,,drop=FALSE]
	if(!is.null(names)) names <- names[DE]
	if(critical.other > critical.primary) warning("critical.other greater than critical.primary")
	otherDE <- (stat > critical.other)
	otherDE[,primary] <- TRUE
	coef[!otherDE] <- NA
	ord <- order(coef[,primary],decreasing=TRUE)

#  Check colours
	if(is.character(low)) low <- col2rgb(low)/255
	if(is.character(high)) high <- col2rgb(high)/255
	col <- rgb( seq(low[1],high[1],len=ncolors), seq(low[2],high[2],len=ncolors), seq(low[3],high[3],len=ncolors) )

#  Truncate coefficients if limit is preset
	if(!is.null(limit))
		if(limit > 0) {
			coef[coef < -limit] <- -limit
			coef[coef > limit] <- limit
		} else
			warning("limit ignored because not positive")

#  Heat plot
	coef <- coef[ord,,drop=FALSE]
	names <- names[ord]
	out <- data.frame(Name=names,coef)
	if(orientation=="portrait") {
		coef <- t(coef)
		coef <- coef[,ng:1,drop=FALSE]
	}
	on.exit(par(old.par))
	if(orientation=="portrait") {
		old.par <- par(mar=cex*c(1,1,4,3))
		image(coef,col=col,xaxt="n",yaxt="n",...)
		cext <- cex*min(1,8/nt)
		mtext(paste(" ",treatments,sep=""),side=3,las=2,at=(cext-1)*0.005+(0:(nt-1))/(nt-1),cex=cext)
		cex <- cex*min(1,40/ng)
		mtext(paste(" ",names,sep=""),side=4,las=2,at=(1-cex)*0.005+((ng-1):0)/(ng-1),cex=cex)
	} else {
		old.par <- par(mar=cex*c(5,6,1,1))
		image(coef,col=col,xaxt="n",yaxt="n",...)
		cext <- cex*min(1,12/nt)
		mtext(paste(treatments," ",sep=""),side=2,las=1,at=(1-cext)*0.005+(0:(nt-1))/(nt-1),cex=cext)
		cex <- cex*min(1,60/ng)
		mtext(paste(names," ",sep=""),side=1,las=2,at=(cex-1)*0.005+(0:(ng-1))/(ng-1),cex=cex)
	}
	invisible(out)
}


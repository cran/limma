#  PLOT DENSITIES

plotDensities<-function(object,log=TRUE,arrays=NULL,singlechannels=NULL,groups=NULL,col=NULL)
#	Plot empirical single-channel densities
#	Original version by Natalie Thorne, 9 September 2003
#	Modified by Gordon Smyth.  Last modified 1 June 2005.  
{
	matDensities<-function(X){
		densXY<-function(Z){
			zd<-density(Z,na.rm=TRUE)
			x<-zd$x
			y<-zd$y
			cbind(x,y)
		}
		out<-apply(X,2,densXY)
		outx<-out[(1:(nrow(out)/2)),]
		outy<-out[(((nrow(out)/2)+1):nrow(out)),]
		list(X=outx,Y=outy)
	}
	
	if(is(object,"MAList")) {
		R <- object$A+object$M/2
		G <- object$A-object$M/2
		if(!log) {
			R <- 2^R
			G <- 2^G
		}
	} else {
		R <- object$R
		G <- object$G
		if(!is.null(object$Rb)) R <- R-object$Rb
		if(!is.null(object$Gb)) G <- G-object$Gb
		if(log) {
			R[R <= 0] <- NA
			G[G <= 0] <- NA
			R <- log(R,2)
			G <- log(G,2)
		}
	}

	if( is.null(arrays) & is.null(singlechannels) ){
		arrays <- 1:(ncol(R))
		x <- cbind(R,G)
		if(is.null(groups)) {
			groups <- c(length(arrays),length(arrays))
			if(is.null(col))
				cols <- rep(c("red","green"),groups)
			if(!is.null(col)) {
				if(length(col)!=2) {
					warning("number of groups=2 not equal to number of col")
					cols<-"black"
				} else {
					cols<-rep(col,groups)
				}
			}
		} else {
			if(!is.null(col)) {
				if(length(as.vector(table(groups)))!=length(col)) {
					warning("number of groups not equal to number of col")
					cols <- col
				} else {
					cols <- col[groups]
				}
			} else {
				warning("warning no cols in col specified for the groups")
				cols <- "black"
			}
		} 
	}else{	
		if(!is.null(singlechannels)) {
			if(!is.null(arrays)) warning("cannot index using arrays AND singlechannels")
			x <- cbind(R,G)[,singlechannels]
			if(is.null(groups)) {
				groups <- c(length(intersect((1:ncol(R)),singlechannels)),
			 length(intersect(((ncol(R)+1):ncol(cbind(G,R))),
									singlechannels)))
				if(is.null(col)) cols <- rep(c("red","green"),groups)
				if(!is.null(col)) {
					if(length(col)!=2) {
						warning("number of groups=2 not equal to number of col")
						cols<-"black"
					} else {
						cols<-rep(col,groups)
					}
				}
			} else {
				if(!is.null(col)) {
					if(length(as.vector(table(groups)))!=length(col)) {
						warning("number of groups not equal to number of col")
						cols <- col
					} else {
						cols <- col[groups]
					}
				} else {
					print("warning no cols in col specified for the groups")
					cols <- "black"
				}
			}
		}else{			
			if(!is.null(arrays)) {
				if(!is.null(singlechannels)) warning("cannot index using arrays AND singlechannels")
				x <- cbind(R[,arrays],G[,arrays])
				if(is.null(groups)) {
					groups <- c(length(arrays),length(arrays))
					if(is.null(col))
						cols <- rep(c("red","green"),groups)
					if(!is.null(col)) {
						if(length(col)!=2) {
							warning("number of groups=2 not equal to number of col")
							cols<-"black"
						} else {
							cols <- rep(col,groups)
						}
					}
				}else{
					if(!is.null(col)) {
						if(length(as.vector(table(groups)))!=length(col)){
							warning("number of groups not equal to number of col")
							cols <- "black"
						}else{
							cols <- col[groups]
						}
					}else{
						warning("warning no cols in col specified for the groups")
						cols <- "black"
					}
				}
			}
		}
	}

	dens.x<-matDensities(x)
	matplot(dens.x$X,dens.x$Y,xlab="Intensity",ylab="Density",main="RG densities",type="l",col=cols,lwd=2,lty=1)
}

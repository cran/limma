#  PLOT DENSITIES

plotDensities<-function(object,log.transform=FALSE,arrays=NULL,singlechannels=NULL,groups=NULL,col=NULL)
#  Plot empirical single-channel densities
#  Natalie Thorne, 9 September 2003
#  Modified by Gordon Smyth, 8 October 2003
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
    
  if(is(object, "MAList")){
    object<- RG.MA(object)
  }

  
  if( ((is.null(arrays))&(is.null(singlechannels))) ){
    arrays <- 1:(ncol(object$R))
    x <- cbind(object$R,object$G)
    if(is.null(groups)){
      groups <- c(length(arrays),length(arrays))
      if(is.null(col))
        colors <- rep(c("red","green"),groups)
      if(!is.null(col)){
        if(length(col)!=2){
          print("number of groups=2 not equal to number of col")
          colors<-"black"
        }else{
          colors<-rep(col,groups)
        }
      }
    }else{
      if(!is.null(col)){
        if(length(as.vector(table(groups)))!=length(col)){
          print("number of groups not equal to number of col")
          colors <- col
        }else{
          colors <- col[groups]
        }                         
      }else{
        print("warning no colors in col specified for the groups")
        colors <- "black"
      }
    } 
  }else{  
    if(!is.null(singlechannels)){
      if(!is.null(arrays))
        print("warning cant index using arrays AND singlechannels")    
      x<-cbind(object$R,object$G)[,singlechannels]
      if(is.null(groups)){
        groups<-c(length(intersect((1:ncol(object$R)),singlechannels)),
       length(intersect(((ncol(object$R)+1):ncol(cbind(object$G,object$R))),
                  singlechannels)))
        if(is.null(col))
          colors <-rep(c("red","green"),groups)
        if(!is.null(col)){
          if(length(col)!=2){
            print("number of groups=2 not equal to number of col")
            colors<-"black"
          }else{
            colors<-rep(col,groups)
          }
        }
      }else{
        if(!is.null(col)){
          if(length(as.vector(table(groups)))!=length(col)){
            print("number of groups not equal to number of col")
            colors <- col
          }else{
            colors <-col[groups]
          }
        }else{
          print("warning no colors in col specified for the groups")
          colors <- "black"
        }
      }
    }else{      
      if(!is.null(arrays)){
        if(!is.null(singlechannels))
          print("warning cant index using arrays AND singlechannels")
        x <- cbind(object$R[, arrays],object$G[,arrays])
        if(is.null(groups)){
          groups <- c(length(arrays),length(arrays))
          if(is.null(col))
            colors <- rep(c("red","green"),groups)
          if(!is.null(col)){
            if(length(col)!=2){
              print("number of groups=2 not equal to number of col")
              colors<-"black"
            }else{
              colors<-rep(col,groups)
            }
          }
        }else{
          if(!is.null(col)){
            if(length(as.vector(table(groups)))!=length(col)){
              print("number of groups not equal to number of col")
              colors<-"black"
            }else{
              colors<-col[groups]
            }
          }else{
            print("warning no colors in col specified for the groups")
            colors <- "black"
          }
        }
      }
    }
  }
  if(log.transform) x <- log(x,2)
  
  dens.x<-matDensities(x)
#  Commented out by GKS  8 Oct 2003
#  XLIM<-c(min(dens.x$X),max(dens.x$X))
#  YLIM<-c(min(dens.x$Y),max(dens.x$Y))
  matplot(dens.x$X,dens.x$Y, xlab = "Intensity", ylab = "Density",
          main = "RG densities",type="l",col=colors,lwd=2,lty=1)
}

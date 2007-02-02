#  limma package /R/mergeScans.R

#  Functions in this file contributed by "Dongseok Choi" <choid@ohsu.edu>
#  Slight modifications for limma version by Gordon Smyth 10 Jan 2007

#-----------------------------------------------------------------------------------------required for mergeScans
.hockey <- function(x,alpha1,beta1,beta2,brk,eps=diff(range(x))/200)
       ## alpha1 is the intercept of the left line segment
       ## beta1 is the slope of the left line segment
       ## beta2 is the slope of the right line segment
       ## brk is location of the break point
       ## 2*eps is the length of the connecting quadratic piece

       ## reference: Bacon & Watts "Estimating the Transition Between
       ## Two Intersecting Straight Lines", Biometrika, 1971

       ## Original function coded by Mary Lindstrom
       ## <lindstro@biostat.wisc.edu> and taken from
       ## S-NEWS Archive (Mon, 24 Apr 2000) available
       ## from (http://lib.stat.cmu.edu/s-news/Burst/15642). 
{
        x1 <- brk-eps
        x2 <- brk+eps
        b <- (x2*beta1-x1*beta2)/(2*eps)
        cc <- (beta2-b)/(2*x2)
        a <- alpha1+beta1*x1-b*x1-cc*x1^2
        alpha2 <- (a + b*x2 + cc*x2^2) - beta2*x2

        lebrk <- (x <= x1)
        gebrk <- (x >= x2)
        eqbrk <- (x > x1 & x < x2)

        result <- rep(0,length(x))
        result[lebrk] <- alpha1 + beta1*x[lebrk]
        result[eqbrk] <- a + b*x[eqbrk] + cc*x[eqbrk]^2
        result[gebrk] <- alpha2 + beta2*x[gebrk]
        result
}

#--------------------------------------------------------------------------------------------required for mergeScans
.mergeScans1 <- function(lowScan, highScan, AboveNoise, outp)
{
         getInitPars <- function(x,y)
         {
           lowLeg <- (x <= quantile(x,probs=0.95))
           init.lsfit <- lm(y~x, subset=lowLeg)$coef
           b0 <- init.lsfit[1]; b1 <- init.lsfit[2];
           b3 <- (y[x==max(x)]-b0)/b1;
           return(c(b0, b1, 0.0, b3))
         }

         findThreshold <- function(x,y)
         {
           n <- length(x)
           pars <- getInitPars(x=x,y=y)
###########################################           
#           if(pars[1]>10)pars[1]<-7
#           if(pars[4]>200)pars[4]<-190
#           print(pars);
######################################
           alpha1 <- pars[1]; beta1 <- pars[2]; beta2 <- pars[3];
           brk <- pars[4];

           fit1 <- nls(y ~ .hockey(x=x,alpha1,beta1,beta2,brk),
                       start=list(alpha1=alpha1,beta1=beta1,beta2=beta2,brk=brk),
                       control=nls.control(maxiter=100))
           resid.scale <- sqrt(fit1$m$deviance()/(n-4))
           yhat <- fit1$m$fitted()
           Coeff <- fit1$m$getPars()
           # cat("\nCoefficients are: ",Coeff,"\n")
           brk <- Coeff[4]; beta2 <- Coeff[3]; beta1 <- Coeff[2]; alpha1 <- Coeff[1];
           Saturated <- y >= -(x - alpha1*beta1 - brk*(beta1*beta1 + 1))/beta1
           # Saturated <- Saturated | (x >= brk)
           high.row <- c(1:length(x))[Saturated]
           High.fit <- alpha1 + beta1*x[high.row]
           return(list(x=x,y=y,yhat=yhat,high.row=high.row,resid.scale = resid.scale,
                       High.fit=High.fit))
         }

         findOutlier <- function(x,y,yhat,outp2)
         {
            n <- length(x)
            xBar <- mean(x)
            Sxx <- (n-1)*var(x)
            resid <- (y-yhat)
            h.ii <- 1/n + (x-xBar)*(x-xBar)/Sxx
            SSE <- sum(resid*resid)
            Rstudent <- resid*sqrt( (n-3)/(SSE*(1-h.ii) - resid^2) )
            pval <- 2*pt(-abs(Rstudent),n-3)
            wksp <- cbind(c(1:n),pval)
            wksp <- wksp[order(wksp[,2]),]
            wksp[1,2] <- min( 1, n*wksp[1,2] )
            i <- 2
            while(i <= n)
            { 
              wksp[i,2] <- min(1, (n+1-i)*wksp[i,2])
              wksp[i,2] <- max(wksp[i,2],wksp[i-1,2])
              i <- i+1
              if (wksp[i-1,2] > outp2) break  
            }
            wksp[,2] <- ((1:n) <= (i-2))*1
            wksp <- wksp[order(wksp[,1]),]
            return(wksp[,2])
         }

    tmp <- findThreshold(x=lowScan,y=highScan)
    high.row <- tmp$high.row
    High.fit <- tmp$High.fit
    yhat <- tmp$yhat
    sigma2 <- tmp$resid.scale^2
    rm(tmp)
    if(outp>0) tmp2 <- findOutlier(x=lowScan[-high.row], y=highScan[-high.row], yhat=yhat[-high.row],outp2=outp)
    else tmp2<-0
    tmp <- matrix(0,nrow=length(yhat),ncol=2)
    tmp[high.row,1] <- 1
    tmp[-high.row,2] <- tmp2

    rm(tmp2)

    dataOut <- cbind(lowScan,highScan,yhat,tmp)
   
    Merged <- yhat*0 - 999
    #################################################################
    ## To Undo a power transformation ( y^lambda) you must use:
    ##     Mean(y) ~ (fit^(1/lambda))*( 1 + sigma^2(1-lambda)/(2*lambda^2*fit^2) )
    ##      var(y) ~ (fit^(2/lambda))*(sigma^2/(lambda^2*fit^2))
    #################################################################
    ok <- dataOut[,4]==1
    Merged[ok] <- High.fit^2 + sigma2
    ok <- (dataOut[,4]==0 & dataOut[,5]==0)
    Merged[ok] <- highScan[ok]^2
    ok <- (dataOut[,4]==0 & dataOut[,5]==1 & AboveNoise==1)
    Merged[ok] <- NA
    ok <- (dataOut[,4]==0 & dataOut[,5]==1 & AboveNoise==0)
    if(sum(ok)!=0) Merged[ok] <- highScan[ok]^2

    rm(yhat,high.row,tmp)

    dataOut <- cbind(dataOut,Merged)
    dimnames(dataOut) <- list(NULL,c("Low","High","yhat","Saturated","Outlier","Merged"))
    return(dataOut)
}

#-------------------------------------------------------------------------------------------handle missing when merging

mergeScansRG <- function(RGlow, RGhigh, AboveNoiseLowG=NULL,AboveNoiseLowR=NULL,outlierp=0.01){
require(limma)
Glow<-RGlow$G
Rlow<-RGlow$R
Ghigh<-RGhigh$G
Rhigh<-RGhigh$R
if (is.null(AboveNoiseLowG)) {AboveNoiseLowG<-matrix(1,nrow=nrow(Glow),ncol=ncol(Glow))}
if (is.null(AboveNoiseLowR)) {AboveNoiseLowR<-matrix(1,nrow=nrow(Rlow),ncol=ncol(Rlow))}

if(ncol(Glow)!=ncol(Ghigh)) {
stop("Number of arrays in low scans and high scans are different")
}
if(ncol(Rlow)!=ncol(Rhigh)) {
stop("Number of arrays in low scans and high scans are different")
}
if(ncol(Ghigh)!=ncol(AboveNoiseLowG)) {
stop("Number of arrays in Green and the number of columns in AboveNoiseHighG are different")
}
if(ncol(Rhigh)!=ncol(AboveNoiseLowR)) {
stop("Number of arrays in Red and the number of columns in AboveNoiseHighR are different")
}


noarrays<-ncol(Glow)
nogenes<-nrow(Glow)
Rmerge<-Gmerge<-matrix(NA, nrow=nogenes,ncol=noarrays)
Rout<-Gout<-matrix(NA, nrow=nogenes,ncol=noarrays)
Rsat<-Gsat<-matrix(NA, nrow=nogenes,ncol=noarrays)
for(i in 1:noarrays){

ok <- apply(sqrt(cbind(Glow[,i],Ghigh[,i])),1,function(x) all(!is.na(x)) )
tmpG <- .mergeScans1(sqrt(Glow[ok,i]),sqrt(Ghigh[ok,i]),AboveNoiseLowG[ok,i],outp=outlierp)
Gmerge[ok,i]<-tmpG[,6]
Gout[ok,i]<-tmpG[,5]
Gsat[ok,i]<-tmpG[,4]

ok <- apply(sqrt(cbind(Rlow[,i],Rhigh[,i])),1,function(x) all(!is.na(x)) )
tmpR <- .mergeScans1(sqrt(Rlow[ok,i]),sqrt(Rhigh[ok,i]),AboveNoiseLowR[ok,i],outp=outlierp)
Rmerge[ok,i]<-tmpR[,6]
Rout[ok,i]<-tmpR[,5]
Rsat[ok,i]<-tmpR[,4]
}

RGmerge<-list(G=Gmerge, R=Rmerge, Gb=RGhigh$Gb, Rb=RGhigh$Gb,other=list(Goutlier=Gout,Routlier=Rout,Gsaturated=Gsat,Rsaturated=Rsat))
new("RGList",RGmerge)

}

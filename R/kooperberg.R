##  KOOPERBERG.R

#  KOOPERBERG BACKGROUND ADUSTMENT FOR GENEPIX DATA

kooperberg <- function(names, fg="mean", bg="median", a=TRUE, layout)
#	Kooperberg Bayesian background correction
#	Matt Ritchie 18 June 2003.
#	Charles Kooperberg contributed 'a' estimation functions (.getas, .varaux1, .varaux2)
#	Modified by Gordon Smyth 16 June 2003.
#	Last modified 7 June 2004.
{
	choices <- c("mean","median")
	fg <- choices[pmatch(fg,choices)]
	bg <- choices[pmatch(bg,choices)]
	meanfg <- (fg=="mean")
	meanbg <- (bg=="mean")
	nslides <- length(names)
	ngenes <- layout$ngrid.c*layout$ngrid.r*layout$nspot.r*layout$nspot.c
	Y <- matrix(0,ngenes,nslides)
	RG <- new("RGList")
	RG$R <- RG$G <- Y
	RG$printer <- layout
	for(i in 1:nslides) {
		temp <- .bayesianAdjustedFG(get(names[i]), meanfg, meanbg, a, layout)
		RG$R[,i] <- temp$R
		RG$G[,i] <- temp$G
	}
	RG
}

.bayesianAdjustedFG <- function(slide, meanfg=TRUE, meanbg=FALSE, a, layout)
#	Matt Ritchie
#	18 June 2003. Last modified 18 September 2004.
{
	ngenes <- dim(slide)[1]
	Y <- rep(0, ngenes)
	RG <- list(R=Y, G=Y)
	if(meanfg) {
		Rfg <- slide[,"F635.Mean"] 
		Gfg <- slide[,"F532.Mean"] 
	} else {
		Rfg <- slide[,"F635.Median"] 
		Gfg <- slide[,"F532.Median"] 
	}
	if(meanbg) {
		Rbg <- slide[,"B635.Mean"] 
		Gbg <- slide[,"B532.Mean"] 
	} else {
		Rbg <- slide[,"B635.Median"] 
		Gbg <- slide[,"B532.Median"] 
	}
	if(a) { 
		aparams <- .getas(slide, layout)
	} else {
		aparams <- c(1,1)
	}
	Rsfg=aparams[2]*slide[,"F635.SD"]/sqrt(slide[,"F.Pixels"]) # column 20
	Rsbg=aparams[2]*slide[,"B635.SD"]/sqrt(slide[,"B.Pixels"]) # column 23
	Gsfg=aparams[1]*slide[,"F532.SD"]/sqrt(slide[,"F.Pixels"]) # column 11
	Gsbg=aparams[1]*slide[,"B532.SD"]/sqrt(slide[,"B.Pixels"]) # column 14
	for(i in 1:ngenes) { 
		if(Rbg[i]>0) {
			RG$R[i] <- .expectedBayesianAdjustedFG(fg = Rfg[i], bg = Rbg[i], sfg = Rsfg[i], sbg = Rsbg[i])	
		} else{
			RG$R[i] <- Rfg[i]
		}
		if(Gbg[i]>0){
	            RG$G[i] <- .expectedBayesianAdjustedFG(fg = Gfg[i], bg = Gbg[i], sfg = Gsfg[i], sbg = Gsbg[i])
		} else {
			RG$G[i] <- Gfg[i]
		}
	}
	RG$R[RG$R>2^16] <- NA
	RG$G[RG$G>2^16] <- NA
	RG
}

.expectedBayesianAdjustedFG <- function(fg, bg, sfg, sbg)
{
	integrate(.numeratorBayesianAdjustedFG, ifelse((fg-bg-4*sqrt(sbg^2+sfg^2))<0, 0, fg-bg-4*sqrt(sbg^2+sfg^2)), 
		ifelse((fg-bg+4*sqrt(sfg^2+sbg^2))<0, 1000, fg-bg+4*sqrt(sfg^2+sbg^2)) , fg=fg, bg=bg, sfg=sfg, sbg=sbg, subdivisions=10000)$value/.denominatorBayesianAdjustedFG(fg, bg, sfg, sbg)
}

.numeratorBayesianAdjustedFG <- function(ut, fg, bg, sfg, sbg)
	ut*exp(dnorm((fg-ut-bg)/sqrt(sfg^2+sbg^2), log=TRUE)+pnorm(((fg-ut)*sbg^2+bg*sfg^2)/(sbg*sfg*sqrt(sfg^2+sbg^2)), log.p=TRUE))

.denominatorBayesianAdjustedFG <- function(fg, bg, sfg, sbg)
{
	sqrt(sfg^2+sbg^2) / sbg * integrate(.normalConvolution,
	ifelse((bg-4*sbg)<0, 0, bg-4*sbg),
	bg+4*sbg, fg=fg, bg=bg, sfg=sfg,
	sbg=sbg, subdivisions=10000)$value
}

.normalConvolution <- function(v, fg, bg, sfg, sbg)
	exp(pnorm((fg-v)/sfg, log.p=TRUE)+dnorm((bg-v)/sbg, log=TRUE))

.getas <- function(x, layout)
{
	b1 <- .varaux1(x,"B532.Mean", layout) # bg G median
	b2 <- .varaux1(x,"B635.Mean", layout) # bg R median
	c1 <- x[, "B532.SD"]/sqrt(x[, "B.Pixels"])
	c2 <- x[, "B635.SD"]/sqrt(x[, "B.Pixels"])
	m1 <- lm(b1 ~ c1 - 1, weights = 1/(c1 + 1))
	m2 <- lm(b2 ~ c2 - 1, weights = 1/(c2 + 1))
	c(m1$coef,m2$coef)
}

# Calculate empirical standard deviation for each spot (based on average of spot and 4 neighbours)
.varaux1 <- function(x, j, layout)
{
	numblocks <- layout$ngrid.c*layout$ngrid.r
	ncols <- layout$nspot.c
	nrows <- layout$nspot.r
	uu <- .varaux2(x, 1, j, ncols, nrows)
	if(numblocks>1) {
		for(i in 2:numblocks) {
		uu <- c(uu, .varaux2(x, i, j, ncols, nrows))
		  }
	}
	uu
}

# Average the standard deviations
.varaux2 <- function(x, i, j, ncols, nrows)
{
	v1 <- x[x[, 1] == i, j]
	v2 <- matrix(v1, ncol=ncols)

# mid grid spot variances
	v4a <- v2[c(-1, -nrows), c(-1, -ncols)]
	v4b <- v2[c(-1, -2), c(-1, -ncols)]
	v4c <- v2[c(-1, -nrows), c(-1, -2)]
	v4d <- v2[c(-(nrows-1), -nrows), c(-1, -ncols)]
	v4e <- v2[c(-1, -nrows), c(-(ncols-1), -ncols)]
	v4x <- cbind(as.vector(v4a), as.vector(v4b), as.vector(v4c),
		 as.vector(v4d), as.vector(v4e))
	VAR <- matrix(0, ncol=ncols, nrow=nrows)
	mid.var <- apply(v4x, 1, FUN=var)
	VAR[2:(nrows-1), 2:(ncols-1)] <- sqrt(mid.var)


# edge spot variances	
# top
	v4a <- v2[1, c(-1, -ncols)]
	v4b <- v2[1, c(-(ncols-1), -ncols)]
	v4c <- v2[2, c(-1, -ncols)]
	v4d <- v2[1, c(-1, -2)]
	v4x <- cbind(as.vector(v4a), as.vector(v4b), as.vector(v4c),
		 as.vector(v4d))
	edge <- apply(v4x, 1, FUN=var)
	VAR[1, 2:(ncols-1)] <- sqrt(edge)

# bottom
	v4a <- v2[nrows, c(-1, -ncols)]
	v4b <- v2[nrows, c(-(ncols-1), -ncols)]
	v4c <- v2[nrows-1, c(-1, -ncols)]
	v4d <- v2[nrows, c(-1, -2)]
	v4x <- cbind(as.vector(v4a), as.vector(v4b), as.vector(v4c),
		 as.vector(v4d))
	
	edge <- apply(v4x, 1, FUN=var)
	VAR[nrows, 2:(ncols-1)] <- sqrt(edge)
	
 # left
	v4a <- v2[c(-1, -nrows), 1]
	v4b <- v2[c(-(nrows-1), -nrows), 1]
	v4c <- v2[c(-1, -nrows), 2]
	v4d <- v2[c(-1, -2), 1]
	v4x <- cbind(as.vector(v4a), as.vector(v4b), as.vector(v4c),
		 as.vector(v4d))
	
	edge <- apply(v4x, 1, FUN=var)
	VAR[2:(nrows-1), 1] <- sqrt(edge)   

 # right
	v4a <- v2[c(-1, -nrows), ncols]
	v4b <- v2[c(-(nrows-1), -nrows), ncols]
	v4c <- v2[c(-1, -nrows), ncols-1]
	v4d <- v2[c(-1, -2), ncols]
	v4x <- cbind(as.vector(v4a), as.vector(v4b), as.vector(v4c),
		 as.vector(v4d))
	
	edge <- apply(v4x, 1, FUN=var)
	VAR[2:(nrows-1), ncols] <- sqrt(edge)   
	
 # corners
	v4x <- cbind(c(v2[1,1], v2[1,ncols], v2[nrows,1], v2[nrows, ncols]),
		  c(v2[1,2], v2[1,ncols-1], v2[nrows,2], v2[nrows, ncols-1]),
		  c(v2[2,1], v2[2,ncols], v2[nrows-1,1], v2[nrows-1, ncols]),
		  c(v2[2,2], v2[2,ncols-1], v2[nrows-1,2], v2[nrows-1, ncols-1]))
	
	corner <- apply(v4x, 1, FUN=var)
	VAR[1,1] <- sqrt(corner[1])
	VAR[1, ncols] <- sqrt(corner[2])
	VAR[nrows, 1] <- sqrt(corner[3])
	VAR[nrows,ncols] <- sqrt(corner[4])
	as.vector(VAR)
}

##  KOOPERBERG.R

#  KOOPERBERG BACKGROUND ADUSTMENT FOR GENEPIX DATA

kooperberg <- function (RG, a = TRUE, layout=RG$printer, verbose=TRUE)
#	Kooperberg Bayesian background correction
#	Matt Ritchie
#	Charles Kooperberg contributed 'a' estimation functions (.getas, .varaux1, .varaux2)
#	Last modified 31 October 2005
{
	if(!is(RG,"RGList")) stop("RG must be an RGList object")
    if(is.null(layout))
	    stop("\nNeed to specify array layout")
    if (is.null(RG$other$"F635 SD") | is.null(RG$other$"B635 SD") | is.null(RG$other$"F532 SD") | is.null(RG$other$"B532 SD") | is.null(RG$other$"B532 Mean") | is.null(RG$other$"B635 Mean") |is.null(RG$other$"F Pixels") | is.null(RG$other$"B Pixels"))
            stop("\nData missing from RG$other: re-run read.maimages with\n other=c(\"F635 SD\",\"B635 SD\",\"F532 SD\",\"B532 SD\",\"B532 Mean\",\"B635 Mean\",\"F Pixels\",\"B Pixels\")")
    nslides <- dim(RG)[2]
    ngenes <- RG$printer$ngrid.c * RG$printer$ngrid.r * RG$printer$nspot.r * RG$printer$nspot.c
    for (i in 1:nslides) {
        temp <- .bayesianAdjustedFG(RG, i, a)
        RG$R[, i] <- temp$R
        RG$G[, i] <- temp$G
	if(verbose)
		{
		cat("Corrected array", i, "\n")
		}
    }
    RG$Rb <- RG$Gb <- NULL 
    RG
}

.bayesianAdjustedFG <- function (RG, k, a = TRUE)
#	Matt Ritchie
#	18 June 2003. Last modified 22 May 2006.
{
    ngenes <- dim(RG)[1] # get number of probes
    Y <- rep(0, ngenes)
    RGmodel <- new("RGList", list(R = Y, G = Y, Rb=NULL, Gb=NULL))

    if (a) {
        aparams <- .getas(RG, k)
    }
    else {
        aparams <- c(1, 1)
    }
    Rsfg = aparams[2] * RG$other$"F635 SD"[,k]/sqrt(RG$other$"F Pixels"[,k])
    Rsbg = aparams[2] * RG$other$"B635 SD"[,k]/sqrt(RG$other$"B Pixels"[,k])
    Gsfg = aparams[1] * RG$other$"F532 SD"[,k]/sqrt(RG$other$"F Pixels"[,k])
    Gsbg = aparams[1] * RG$other$"B532 SD"[,k]/sqrt(RG$other$"B Pixels"[,k])
    for (i in 1:ngenes) {
        if (RG$R[i,k] > 0 & Rsbg[i] > 0) {
            RGmodel$R[i] <- .expectedBayesianAdjustedFG(fg = RG$R[i,k],
                bg = RG$Rb[i,k], sfg = Rsfg[i], sbg = Rsbg[i])
        }
        else {
            RGmodel$R[i] <- RG$R[i,k]
        }
        if (RG$G[i,k] > 0 & Gsbg[i] > 0) {
            RGmodel$G[i] <- .expectedBayesianAdjustedFG(fg = RG$G[i,k],
                bg = RG$Gb[i,k], sfg = Gsfg[i], sbg = Gsbg[i])
        }
        else {
            RGmodel$G[i] <- RG$G[i,k]
        }
    }
    RGmodel$R[RGmodel$R > 2^16] <- NA
    RGmodel$G[RGmodel$G > 2^16] <- NA
    RGmodel
}

.getas <- function (RG, j)
{
    b1 <- .varaux1(RG$other$"B532 Mean"[,j], RG$printer)
    b2 <- .varaux1(RG$other$"B635 Mean"[,j], RG$printer)
    c1 <- RG$other$"B532 SD"[,j]/sqrt(RG$other$"B Pixels"[,j])
    c2 <- RG$other$"B635 SD"[,j]/sqrt(RG$other$"B Pixels"[,j])
    m1 <- lm(b1 ~ c1 - 1, weights = 1/(c1 + 1))
    m2 <- lm(b2 ~ c2 - 1, weights = 1/(c2 + 1))
    c(m1$coef, m2$coef)
}

# Calculate empirical standard deviation for each spot (based on average of spot and 4 neighbours)
.varaux1 <- function (bg, layout)
{
    numblocks <- layout$ngrid.c * layout$ngrid.r
    block <- rep(1:numblocks, each=layout$nspot.r*layout$nspot.c)
    uu <- .varaux2(bg, block, 1, layout$nspot.c, layout$nspot.r)
    if (numblocks > 1) {
        for (i in 2:numblocks) {
            uu <- c(uu, .varaux2(bg, block, i, layout$nspot.c, layout$nspot.r))
        }
    }
    uu
}

# Average the standard deviations
.varaux2 <- function (bg, block, i, ncols, nrows)
{
    v1 <- bg[block == i]
    v2 <- matrix(v1, ncol = ncols)

# mid grid spot variances    
    v4a <- v2[c(-1, -nrows), c(-1, -ncols)]
    v4b <- v2[c(-1, -2), c(-1, -ncols)]
    v4c <- v2[c(-1, -nrows), c(-1, -2)]
    v4d <- v2[c(-(nrows - 1), -nrows), c(-1, -ncols)]
    v4e <- v2[c(-1, -nrows), c(-(ncols - 1), -ncols)]
    v4x <- cbind(as.vector(v4a), as.vector(v4b), as.vector(v4c),
        as.vector(v4d), as.vector(v4e))
    VAR <- matrix(0, ncol = ncols, nrow = nrows)
    mid.var <- apply(v4x, 1, FUN = var)
    VAR[2:(nrows - 1), 2:(ncols - 1)] <- sqrt(mid.var)

# edge spot variances	
# top
    v4a <- v2[1, c(-1, -ncols)]
    v4b <- v2[1, c(-(ncols - 1), -ncols)]
    v4c <- v2[2, c(-1, -ncols)]
    v4d <- v2[1, c(-1, -2)]
    v4x <- cbind(as.vector(v4a), as.vector(v4b), as.vector(v4c),
        as.vector(v4d))
    edge <- apply(v4x, 1, FUN = var)
    VAR[1, 2:(ncols - 1)] <- sqrt(edge)

# bottom
    v4a <- v2[nrows, c(-1, -ncols)]
    v4b <- v2[nrows, c(-(ncols - 1), -ncols)]
    v4c <- v2[nrows - 1, c(-1, -ncols)]
    v4d <- v2[nrows, c(-1, -2)]
    v4x <- cbind(as.vector(v4a), as.vector(v4b), as.vector(v4c),
        as.vector(v4d))
    edge <- apply(v4x, 1, FUN = var)
    VAR[nrows, 2:(ncols - 1)] <- sqrt(edge)

# left
    v4a <- v2[c(-1, -nrows), 1]
    v4b <- v2[c(-(nrows - 1), -nrows), 1]
    v4c <- v2[c(-1, -nrows), 2]
    v4d <- v2[c(-1, -2), 1]
    v4x <- cbind(as.vector(v4a), as.vector(v4b), as.vector(v4c),
        as.vector(v4d))
    edge <- apply(v4x, 1, FUN = var)
    VAR[2:(nrows - 1), 1] <- sqrt(edge)

# right
    v4a <- v2[c(-1, -nrows), ncols]
    v4b <- v2[c(-(nrows - 1), -nrows), ncols]
    v4c <- v2[c(-1, -nrows), ncols - 1]
    v4d <- v2[c(-1, -2), ncols]
    v4x <- cbind(as.vector(v4a), as.vector(v4b), as.vector(v4c),
        as.vector(v4d))
    edge <- apply(v4x, 1, FUN = var)
    VAR[2:(nrows - 1), ncols] <- sqrt(edge)

# corners   
    v4x <- cbind(c(v2[1, 1], v2[1, ncols], v2[nrows, 1], v2[nrows,
        ncols]), c(v2[1, 2], v2[1, ncols - 1], v2[nrows, 2],
        v2[nrows, ncols - 1]), c(v2[2, 1], v2[2, ncols], v2[nrows -
        1, 1], v2[nrows - 1, ncols]), c(v2[2, 2], v2[2, ncols -
        1], v2[nrows - 1, 2], v2[nrows - 1, ncols - 1]))
    corner <- apply(v4x, 1, FUN = var)
    VAR[1, 1] <- sqrt(corner[1])
    VAR[1, ncols] <- sqrt(corner[2])
    VAR[nrows, 1] <- sqrt(corner[3])
    VAR[nrows, ncols] <- sqrt(corner[4])
    as.vector(VAR)
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

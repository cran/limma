#  SUBSET AND COMBINE DATA SETS

assign("[.RGList",
function(object, i, j, ...) {
#  Subsetting for RGList objects
#  Gordon Smyth
#  29 June 2003.  Last modified 22 December 2005.

	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	oc <- names(object$other)
	if(missing(i))
		if(missing(j))
			return(object)
		else {
			object$R <- object$R[,j,drop=FALSE]
			object$G <- object$G[,j,drop=FALSE]
			object$Rb <- object$Rb[,j,drop=FALSE]
			object$Gb <- object$Gb[,j,drop=FALSE]
			object$weights <- object$weights[,j,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			for(k in oc) object$other[[k]] <- object$other[[k]][,j,drop=FALSE]
		}
	else {
		if(is.character(i)) {
			i <- match(i,rownames(object))
			i <- i[!is.na(i)]
		}
		if(missing(j)) {
			object$R <- object$R[i,,drop=FALSE]
			object$G <- object$G[i,,drop=FALSE]
			object$Rb <- object$Rb[i,,drop=FALSE]
			object$Gb <- object$Gb[i,,drop=FALSE]
			object$weights <- object$weights[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			for(k in oc) object$other[[k]] <- object$other[[k]][i,,drop=FALSE]
		} else {
			object$R <- object$R[i,j,drop=FALSE]
			object$G <- object$G[i,j,drop=FALSE]
			object$Rb <- object$Rb[i,j,drop=FALSE]
			object$Gb <- object$Gb[i,j,drop=FALSE]
			object$weights <- object$weights[i,j,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			for(k in oc) object$other[[k]] <- object$other[[k]][i,j,drop=FALSE]
		}
	}
	object
})

assign("[.MAList",
function(object, i, j, ...) {
#  Subsetting for MAList objects
#  Gordon Smyth
#  29 June 2003.  Last modified 22 Dec 2005.

	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	other <- names(object$other)
	if(missing(i))
		if(missing(j))
			return(object)
		else {
			object$M <- object$M[,j,drop=FALSE]
			object$A <- object$A[,j,drop=FALSE]
			object$weights <- object$weights[,j,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			if(!is.null(object$design)) {
				object$design <- as.matrix(object$design)[j,,drop=FALSE]
				if(!is.fullrank(object$design)) warning("subsetted design matrix is singular",call.=FALSE)
			}
			for(a in other) object$other[[a]] <- object$other[[a]][,j,drop=FALSE]
		}
	else {
		if(is.character(i)) {
			i <- match(i,rownames(object))
			i <- i[!is.na(i)]
		}
		if(missing(j)) {
			object$M <- object$M[i,,drop=FALSE]
			object$A <- object$A[i,,drop=FALSE]
			object$weights <- object$weights[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			for(a in other) object$other[[a]] <- object$other[[a]][i,,drop=FALSE]
		} else {
			object$M <- object$M[i,j,drop=FALSE]
			object$A <- object$A[i,j,drop=FALSE]
			object$weights <- object$weights[i,j,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$targets <- object$targets[j,,drop=FALSE]
			if(!is.null(object$design)) {
				object$design <- as.matrix(object$design)[j,,drop=FALSE]
				if(!is.fullrank(object$design)) warning("subsetted design matrix is singular",call.=FALSE)
			}
			for(a in other) object$other[[a]] <- object$other[[a]][i,j,drop=FALSE]
		}
	}
	object
})

assign("[.MArrayLM",
function(object, i, j, ...)
#  Subsetting for MArrayLM objects
#  Gordon Smyth
#  26 April 2005. Last modified 30 August 2006.
{
	if(nargs() != 3) stop("Two subscripts required",call.=FALSE)
	if(!is.null(object$coefficients)) object$coefficients <- as.matrix(object$coefficients)
	if(!is.null(object$stdev.unscaled)) object$stdev.unscaled <- as.matrix(object$stdev.unscaled)
	if(!is.null(object$weights)) object$weights <- as.matrix(object$weights)
	if(!is.null(object$p.value)) object$p.value <- as.matrix(object$p.value)
	if(!is.null(object$lods)) object$lods <- as.matrix(object$lods)
	if(!is.null(object$targets)) object$targets <- as.data.frame(object$targets)
	if(!is.null(object$cov.coefficients)) object$cov.coefficients <- as.matrix(object$cov.coefficients)
	if(!is.null(object$contrasts)) object$contrasts <- as.matrix(object$contrasts)
	if(is.null(object$contrasts) && !is.null(object$coefficients)) {
		object$contrasts <- diag(ncol(object$coefficients))
		rownames(object$contrasts) <- colnames(object$contrasts) <- colnames(object$coefficients)
	}
	if(!is.null(object$genes)) object$genes <- as.data.frame(object$genes)
	if(missing(i)) {
		if(missing(j))
			return(object)
		else {
			object$coefficients <- object$coefficients[,j,drop=FALSE]
			object$stdev.unscaled <- object$stdev.unscaled[,j,drop=FALSE]
			object$t <- object$t[,j,drop=FALSE]
			object$weights <- object$weights[,j,drop=FALSE]
			object$p.value <- object$p.value[,j,drop=FALSE]
			object$lods <- object$lods[,j,drop=FALSE]
			object$cov.coefficients <- object$cov.coefficients[j,j,drop=FALSE]
			object$contrasts <- object$contrasts[,j,drop=FALSE]
			object$var.prior <- object$var.prior[j]
		}
	} else {
		if(is.character(i)) {
			i <- match(i,rownames(object))
			i <- i[!is.na(i)]
		}
		if(missing(j)) {
			object$coefficients <- object$coefficients[i,,drop=FALSE]
			object$stdev.unscaled <- object$stdev.unscaled[i,,drop=FALSE]
			object$t <- object$t[i,,drop=FALSE]
			object$weights <- object$weights[i,,drop=FALSE]
			object$p.value <- object$p.value[i,,drop=FALSE]
			object$lods <- object$lods[i,,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
		} else {
			object$coefficients <- object$coefficients[i,j,drop=FALSE]
			object$stdev.unscaled <- object$stdev.unscaled[i,j,drop=FALSE]
			object$t <- object$t[i,j,drop=FALSE]
			object$weights <- object$weights[i,j,drop=FALSE]
			object$p.value <- object$p.value[i,j,drop=FALSE]
			object$lods <- object$lods[i,j,drop=FALSE]
			object$genes <- object$genes[i,,drop=FALSE]
			object$cov.coefficients <- object$cov.coefficients[j,j,drop=FALSE]
			object$contrasts <- object$contrasts[,j,drop=FALSE]
			object$var.prior <- object$var.prior[j]
		}
		object$df.residual <- object$df.residual[i]
		object$sigma <- object$sigma[i]
		object$s2.post <- object$s2.post[i]
		object$Amean <- object$Amean[i]
		object$F <- object$F[i]
		object$F.p.value <- object$F.p.value[i]
	}
	object
})

cbind.RGList <- function(..., deparse.level=1) {
#  Combine RGList objects assuming same genelists
#  Gordon Smyth
#  27 June 2003. Last modified 6 Nov 2005.

	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	other <- names(objects[[1]]$other)
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$R <- cbind(out$R,objects[[i]]$R)
		out$G <- cbind(out$G,objects[[i]]$G)
		out$Rb <- cbind(out$Rb,objects[[i]]$Rb)
		out$Gb <- cbind(out$Gb,objects[[i]]$Gb)
		out$weights <- cbind(out$weights,objects[[i]]$weights)
		out$targets <- rbind(out$targets,objects[[i]]$targets)
		for (a in other) out$other[[a]] <- cbind(out$other[[a]],objects[[i]]$other[[a]])
	}
	out
}

cbind.MAList <- function(..., deparse.level=1) {
#  Combine MAList objects assuming same genelists
#  Gordon Smyth
#  27 June 2003. Last modified 6 Nov 2005.

	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	other <- names(objects[[1]]$other)
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$M <- cbind(out$M,objects[[i]]$M)
		out$A <- cbind(out$A,objects[[i]]$A)
		out$weights <- cbind(out$weights,objects[[i]]$weights)
		out$targets <- rbind(out$targets,objects[[i]]$targets)
		for (a in other) out$other[[a]] <- cbind(out$other[[a]],objects[[i]]$other[[a]])
	}
	out
}

rbind.RGList <- function(..., deparse.level=1) {
#  Combine RGList objects assuming same array lists
#  Gordon Smyth
#  6 Dec 2003. Last modified 6 Nov 2005.

	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	other <- names(objects[[1]]$other)
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$R <- rbind(out$R,objects[[i]]$R)
		out$G <- rbind(out$G,objects[[i]]$G)
		out$Rb <- rbind(out$Rb,objects[[i]]$Rb)
		out$Gb <- rbind(out$Gb,objects[[i]]$Gb)
		out$weights <- rbind(out$weights,objects[[i]]$weights)
		out$genes <- rbind(out$genes,objects[[i]]$genes)
		for (a in other) out$other[[a]] <- rbind(out$other[[a]],objects[[i]]$other[[a]])
	}
	out
}

rbind.MAList <- function(..., deparse.level=1) {
#  Combine MAList objects assuming same array lists
#  Gordon Smyth
#  7 Dec 2003. Last modified 6 Nov 2005.

	objects <- list(...)
	nobjects <- length(objects)
	out <- objects[[1]]
	other <- names(objects[[1]]$other)
	if(nobjects > 1)
	for (i in 2:nobjects) {
		out$M <- rbind(out$M,objects[[i]]$M)
		out$A <- rbind(out$A,objects[[i]]$A)
		out$weights <- rbind(out$weights,objects[[i]]$weights)
		out$genes <- rbind(out$genes,objects[[i]]$genes)
		for (a in other) out$other[[a]] <- rbind(out$other[[a]],objects[[i]]$other[[a]])
	}
	out
}

makeUnique <- function(x) {
#  Add characters to the elements of a character vector to make all values unique
#  Gordon Smyth
#  10 April 2003

	x <- as.character(x)
	tab <- table(x)
	tab <- tab[tab>1]
	lentab <- length(tab)
	if(lentab > 0) {
		u <- names(tab)
		for (i in 1:lentab) {
			n <- tab[i]
			x[x==u[i]] <- paste(x[x==u[i]],formatC(1:n,width=1+floor(log(n,10)),flag="0"),sep="")
		}
	}
	x
}

merge.RGList <- function(x,y,...) {
#  Merge RGList y into x aligning by row names
#  Gordon Smyth
#  11 April 2003.  Last modified 28 Oct 2005.

	if(!is(y,"RGList")) stop("both x and y must be RGList objects")
	genes1 <- rownames(x$R)
	if(is.null(genes1)) genes1 <- rownames(x$G)
	if(is.null(genes1)) genes1 <- x$genes$ID
	genes2 <- rownames(y$R)
	if(is.null(genes2)) genes2 <- rownames(y$G)
	if(is.null(genes2)) genes2 <- y$genes$ID
	if(is.null(genes1) || is.null(genes2)) stop("Need row names to align on") 

	fields1 <- names(x)
	fields2 <- names(y)
	if(!identical(fields1,fields2)) stop("The two RGLists have different components")

	ord2 <- match(makeUnique(genes1), makeUnique(genes2))
	cbind(x,y[ord2,])
}

merge.MAList <- function(x,y,...) {
#  Merge MAList y into x aligning by row names
#  Gordon Smyth
#  7 May 2004.  Last modified 28 Oct 2005.

	if(!is(y,"MAList")) stop("both x and y must be MAList objects")
	genes1 <- rownames(x$M)
	if(is.null(genes1)) genes1 <- rownames(x$A)
	if(is.null(genes1)) genes1 <- x$genes$ID
	genes2 <- rownames(y$M)
	if(is.null(genes2)) genes2 <- rownames(y$A)
	if(is.null(genes2)) genes2 <- y$genes$ID
	if(is.null(genes1) || is.null(genes2)) stop("Need row names to align on") 

	fields1 <- names(x)
	fields2 <- names(y)
	if(!identical(fields1,fields2)) stop("The two MALists have different components")

	ord2 <- match(makeUnique(genes1), makeUnique(genes2))
	cbind(x,y[ord2,])
}


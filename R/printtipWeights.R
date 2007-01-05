printtipWeights <- function(object, design = NULL, weights = NULL, method="genebygene", layout = object$printer, maxiter=50, tol = 1e-10, trace = FALSE)
#	Compute print-tip quality weights
#	Matt Ritchie
#	21 Sep 2006. Last revised 4 Dec 2006.
{
	M <- NULL
	if (is(object, "MAList") || is(object, "list")) {
		M <- object$M
		if (missing(design) && !is.null(object$design))
			design <- object$design
		if (missing(weights) && !is.null(object$weights))
			weights <- object$weights
	} else {
		if (is(object, "marrayNorm")) {
			M <- object@maM
			if (missing(weights) && length(object@maW))
				weights <- object@maW
		} else {
			if (is(object, "exprSet"))
				M <- object@exprs
		}
	}

	if (is.null(M))
		M <- as.matrix(object)
	if (is.null(design))
		design <- matrix(1, ncol(M), 1)
	design <- as.matrix(design)

	if(mode(design) != "numeric") stop("design must be a numeric matrix")
	    ne <- nonEstimable(design)
	if(!is.null(ne))
            cat("Coefficients not estimable:",paste(ne,collapse=" "),"\n")
        p <- ncol(design)
#        cols <- seq(1:p)
        QR <- qr(design)
        nparams <- QR$rank # ncol(design)
        nprobes <- nrow(M)
	narrays <- ncol(M)
        if(is.null(layout)) stop("printer layout information must be specified")
        ngr <- layout$ngrid.r
	ngc <- layout$ngrid.c
	nspots <- layout$nspot.r * layout$nspot.c
	nprobes2 <- ngr*ngc*nspots
	if(nprobes2 != nprobes) stop("printer layout information does not match M row dimension")
        printtipwts <- matrix(0, ngr*ngc, narrays)
	# Set up design matrix for array variance model
	Z <- contr.sum(narrays)

	spots <- 1:nspots
	gridnum <- 1
	for (gridr in 1:ngr) {
		for (gridc in 1:ngc) {	
			# Intialise sub-array variances to zero
			arraygammas <- rep(0, (narrays-1))

			method <- match.arg(method,c("genebygene","reml"))
			switch(method, genebygene = {  # Estimate sub-array variances via gene-by-gene update algorithm
				Zinfo <- 10*(narrays-nparams)/narrays*crossprod(Z, Z)
				for(i in spots) {
                                        if(!all(is.finite(arraygammas)))
                                                stop("convergence problem at gene ", i, " in print-tip grid row ", gridr, ", column ", gridc,  ": array weights not estimable")

				 	vary <- exp(Z%*%arraygammas)

					if(!is.null(weights)) {  # combine spot weights with running weights 
						if(max(weights[i,], na.rm=TRUE) > 1) {
							weights[i,] <- weights[i,]/max(weights[i,])
						}
						w <- as.vector(1/vary*weights[i,])
					} else {
						w <- as.vector(1/vary)
					}
					y <- as.vector(M[i,])
					obs <- is.finite(y) & w!=0
					if (sum(obs) > 1) {
						if(sum(obs) == narrays)	{
							X <- design
						} else {  # remove missing/infinite values
							X <- design[obs, , drop = FALSE]
							y <- y[obs]
							vary <- vary[obs]
							Z2 <- Z[obs,]
						}
		
						out <- lm.wfit(X, y, w[obs])
						d <- rep(0, narrays)
						d[obs] <- w[obs]*out$residuals^2
 						s2 <- sum(d[obs])/out$df.residual
                                                Q <- qr.Q(out$qr)
                                                if(ncol(Q)!=out$rank)
                                                        Q <- Q[,-((out$rank+1):ncol(Q)),drop=FALSE]
						h <- rep(1, narrays)
						h[obs] <- rowSums(Q^2)
						Agam <- crossprod(Z, (1-h)*Z)
						Agam.del <- crossprod(t(rep(h[narrays], narrays-1)-h[1:(length(narrays)-1)]))
						Agene.gam <- (Agam - 1/out$df.residual*Agam.del)
						if(is.finite(sum(Agene.gam)) && sum(obs) == narrays) {
							Zinfo <- Zinfo + Agene.gam
							R <- chol(Zinfo)
							Zinfoinv <- chol2inv(R)
							zd <- d/s2 - 1 + h
							Zzd <- crossprod(Z, zd)
							gammas.iter <- Zinfoinv%*%Zzd
							arraygammas <- arraygammas + gammas.iter
						}
						if(is.finite(sum(Agene.gam)) && sum(obs) < narrays && sum(obs) > 2) { 
							Zinfo <- Zinfo + Agene.gam
							A1 <- (diag(1, sum(obs))-1/sum(obs)*matrix(1, sum(obs), sum(obs)))%*%Z2
							A1 <- A1[-sum(obs),] # remove last row
							R <- chol(Zinfo)
							Zinfoinv <- chol2inv(R)
							zd <- d/s2 - 1 + h
							Zzd <- A1%*%crossprod(Z, zd)
							Zinfoinv.A1 <- A1%*%Zinfoinv%*%t(A1)
							alphas.old <- A1%*%arraygammas
							alphas.iter <- Zinfoinv.A1%*%Zzd
							alphas.new <- alphas.old + alphas.iter
							Us <- rbind(diag(1, sum(obs)-1), -1)
							Usalphas <- Us%*%(alphas.new-alphas.old)
							Usgammas <- Z%*%arraygammas
							Usgammas[obs] <- Usgammas[obs] + Usalphas
							arraygammas <- Usgammas[1:(narrays-1)]
						}
					}
				}
			}, reml = {  # Estimate sub-array variances via reml
#				const <- narrays * log(2 * pi)
				iter <- 0
				dev <- 0
				repeat {
#			devold <- dev
#			dev <- 0
				iter <- iter + 1
				zd <- matrix(0, narrays, 1)
				sum1minush <- matrix(0, narrays, 1)
				K <- matrix(0, nprobes, narrays)
	
				for(i in spots) {
					vary <- exp(Z%*%arraygammas)
	
					if(!is.null(weights)) {  # combine spot weights with running weights 
						if(max(weights[i,], na.rm=TRUE) > 1) {
							weights[i,] <- weights[i,]/max(weights[i,])
						}
						w <- as.vector(1/vary*weights[i,])
					} else {
						w <- as.vector(1/vary)
					}
	
					y <- as.vector(M[i,])
					obs <- is.finite(y) & w!=0
					n <- sum(obs)
					if (n > 0) {
						if(n == narrays)	{
							X <- design
							#Z2 <- Z
						} else {  # remove missing/infinite values
							X <- design[obs, , drop = FALSE]
							y <- y[obs]
							vary <- vary[obs]
							w <- w[obs]
							const <- sum(obs) * log(2 * pi)
						}
						# cat(i)
						out <- lm.wfit(X, y, w)
						d <- rep(0, narrays)
						d[obs] <- w*out$residuals^2
						s2 <- sum(d[obs])/out$df.residual
                                                Q <- qr.Q(out$qr)
                                                if(ncol(Q)!=out$rank)
                                                         Q <- Q[,-((out$rank+1):ncol(Q)),drop=FALSE]
						h <- rowSums(Q^2)
						zd[obs] <- zd[obs] + d[obs]/s2 - 1 + h
						sum1minush[obs,1] <- sum1minush[obs,1] + 1-h
						K[i,][obs] <- as.vector(h[n]-h)
#						dev <- dev + sum(d[obs]/vary) + sum(log(vary)) + const + 2 * log(prod(abs(diag(out$qr$qr))))
					}
				}
				Zzd <- crossprod(Z, zd)
				Zinfo <- diag(sum1minush[1:(narrays-1)]) + sum1minush[narrays] - crossprod(K[,-narrays])/out$df.residual # (narrays-nparams)
				R <- chol(Zinfo)
				Zinfoinv <- chol2inv(R)
				gammas.iter <- Zinfoinv%*%Zzd
				arraygammas <- arraygammas + gammas.iter
#				arrayw <- drop(exp(Z %*% (-arraygammas)))
				x2 <- crossprod(Zzd, gammas.iter) / narrays
	
				if(trace)
					cat("Iter =", iter, " X2 =", x2, " Print-tip grid row ", gridr, 
							" column ", gridc, " gammas: ", arraygammas, "\n")

                                if(!all(is.finite(arraygammas)))
                                        stop("convergence problem at iteration ", iter, ": array weights not estimable")
	
#				if (dev < devold - 1e-50)
#					break

				if (x2  < tol)
					break

				if (iter == maxiter)	{
					warning("Maximum iterations ", maxiter, " reached", sep="")
					break
				}
				}
	
			})
#			matrix(rep(1/exp(Z%*%arraygammas), each=nspots), nprobes, narrays)
			printtipwts[gridnum, ] <- drop(exp(Z %*% (-arraygammas)))
			spots <- spots + nspots
#			arrays <- arrays + narrays
			gridnum <- gridnum + 1
#			cat("\n")
		}
	}
	# Make a matrix of print-tip weights for use in lmFit()
	ngrids <- ngr*ngc
	wts <- matrix(0, nspots*ngrids, narrays)	
	spots <- 1:nspots
	for(i in 1:ngrids) {
		for(j in 1:narrays) {
			wts[spots,j] <- printtipwts[i,j] 
			}
	spots <- spots + nspots
	}
	wts
}

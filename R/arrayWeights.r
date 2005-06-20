arrayWeights <- function(object, design = NULL, weights = NULL)
#	Compute array quality weights
#	Matt Ritchie
#	7 Feb 2005. Last revised 5 May 2005.
{
	M <- NULL
	if (is(object, "MAList") || is(object, "list")) {
	        M <- object$M
		if (missing(design) && !is.null(object$design)) 
		            design <- object$design
		if (missing(weights) && !is.null(object$weights)) 
		            weights <- object$weights
		}
	else {
	        if (is(object, "marrayNorm")) {
	            M <- object@maM
	            if (missing(weights) && length(object@maW)) 
	                weights <- object@maW
		        }
   	            else {
		          if (is(object, "exprSet")) {
                		M <- object@exprs
            			}
            	    else {
            		  if (is(object, "PLMset")) {
                  	  M <- object@chip.coefs
                  	  if (missing(weights) && length(object@se.chip.coefs)) 
                  	  	weights <- 1/pmax(object@se.chip.coefs, 1e-05)^2
                		}
            		  }
        	    }
    		}
    	if (is.null(M)) 
    	    M <- as.matrix(object)
    	if (is.null(design)) 
    	    design <- matrix(1, ncol(M), 1)
    	design <- as.matrix(design)

	params <- ncol(design)
	ngenes <- dim(M)[1]
	narrays <- dim(M)[2]
	
	# Set up design matrix for array variance model
	arrays <- as.factor(1:narrays)
	contrasts(arrays) <- contr.sum(levels(arrays))
	Z <- model.matrix(~as.factor(arrays))
	Z <- Z[,-1]

	# Intialise array variances to zero
	arraygammas <- rep(0, (narrays-1))
	
	# Estimate array variances via one-step update
	Zinfo <- 10*(1-1/narrays)*crossprod(Z, Z)
	for(i in 1:ngenes) {
	 	vary <- exp(Z%*%arraygammas)
		
		if(!is.null(weights)) {		# combine spot weights with running weights 
			if(max(weights[i,], na.rm=TRUE) > 1) {
				weights[i,] <- weights[i,]/max(weights[i,])
				}
			w <- as.vector(1/vary*weights[i,])
			}
		else {
			w <- as.vector(1/vary)
			}
		y <- as.vector(M[i,])
		obs <- is.finite(y) & w!=0
		if (sum(obs) > 1) {		
			if(sum(obs) == narrays)	{
				X <- design
				}
			else {		# remove missing/infinite values
				X <- design[obs, , drop = FALSE]
				y <- y[obs]
				vary <- vary[obs]
				Z2 <- Z[obs,]
				}
							
			out <- lm.wfit(X, y, w[obs])
			d <- rep(0, narrays)
			d[obs] <- w[obs]*out$residuals^2
			Q <- qr.Q(out$qr)
			nparams <- dim(Q)[2]
			h <- rep(1, narrays)
			h[obs] <- as.vector(Q^2 %*% array(1, c(nparams, 1)))
			Agam <- crossprod(Z, (1-h)*Z)
			Agam.del <- crossprod(t(rep(h[narrays], narrays-1)-h[1:(length(narrays)-1)]))
			Agene.gam <- (Agam - 1/(narrays-params)*Agam.del)
			if(is.finite(sum(Agene.gam)) && sum(obs) == narrays) {
				Zinfo <- Zinfo + Agene.gam
				R <- chol(Zinfo)
				Zinfoinv <- chol2inv(R)
				zd <- d - 1 + h
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
				zd <- d - 1 + h
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
	matrix(rep(1/exp(Z%*%arraygammas), each=ngenes), ngenes, narrays)
}

arrayWeights <- function(object, design = NULL, weights = NULL)
#	Compute array quality weights
#	Matt Ritchie
#	7 Feb 2005.
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
	Zinfo <- 10*crossprod(Z)
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
		if (sum(obs) > 0) {		
			if(sum(obs) == narrays)	{
				X <- design
				Z2 <- Z
				}
			else {		# remove missing/infinite values
				X <- design[obs, , drop = FALSE]
				y <- y[obs]
				vary <- vary[obs]
				w <- w[obs]
				Z2 <- Z[obs, , drop = FALSE] 
				}
							
			out <- lm.wfit(X, y, w)
			d <- w*out$residuals^2
			Q <- qr.Q(out$qr)
			nparams <- dim(Q)[2]
			h <- as.vector(Q^2 %*% array(1, c(nparams, 1)))
			Zinfo <- Zinfo + crossprod(Z2, (1-h)*Z2)
			R <- chol(Zinfo)
			zd <- d - 1 + h
			Zzd <- crossprod(Z2, zd)
			
			gammas.iter <- backsolve(R, backsolve(R, Zzd, transpose = TRUE))
			arraygammas <- arraygammas + gammas.iter
			}
		}
	matrix(rep(1/exp(Z%*%arraygammas), each=ngenes), ngenes, narrays)
}

fit.normexp <- function(foreground,background)
#	Fit background=normal + signal=exponential model using BFGS.
#	Jeremy Silver.
#	21 Jan 2005. Last modified 8 Feb 2005.
{
	#the following are rough estimates (starting values) for parameters beta, alpha and sigma.
	mu <- quantile(foreground,0.03, na.rm = TRUE, names = FALSE)
	beta <- mu - mean(background, na.rm = TRUE)
	sigma <- sqrt(mean((foreground[foreground<mu]-mu)^2, na.rm = TRUE))
	alpha <- mean(foreground,na.rm = TRUE) - mu
	z<-0												# This is an indicator that tells us if we need to rescale or not
		
	if (is.infinite(sumloglik(c(beta,log(sigma),log(alpha)),foreground,background))){	#this step prevents calculating infinite sums of loglikelihoods
		foreground<-foreground/1000						# Make the values a bit smaller
		background<-background/1000	
		z<-1											# remind ourselves to rescale later
		mu <- quantile(foreground,0.03, na.rm = TRUE, names = FALSE)
		beta <- mu - mean(background, na.rm = TRUE)
		sigma <- sqrt(mean((foreground[foreground<mu]-mu)^2, na.rm = TRUE))
		alpha <- mean(foreground,na.rm = TRUE) - mu
		if(is.infinite(sumloglik(c(beta,log(alpha),log(sigma)),foreground,background))){
			beta<-1										# If the new estimates are still infinite we give some
			alpha<-exp(1)								# uneducated guess and hope that optim will manage
			sigma<-exp(1)
		}
	}
	
	Results <- optim(c(beta,log(sigma),log(alpha)), sumloglik, gr=grsumloglik, method=c("BFGS"), lower = -Inf, upper = Inf, control = list(), hessian = FALSE, foreground, background)
	
	if (z==0){											# the case where we don't need to rescale
		out <- list(beta=Results[[1]][1], sigma=exp(Results[[1]][2]), alpha=exp(Results[[1]][3]),m2loglik=Results[[2]], convergence=Results[[4]])
	} else {											# in which we have to rescale
		out <- list(beta=1000*Results[[1]][1], sigma=1000*exp(Results[[1]][2]), alpha=1000*exp(Results[[1]][3]),m2loglik=m2loglik.normexp(c(1000*Results[[1]][1],log(1000*exp(Results[[1]][2])),log(1000*exp(Results[[1]][3]))),1000*foreground,1000*background),convergence=Results[[4]])
	}
	out	
}

grsumloglik <- function(theta,foreground,background)
#	Gradient of norm-exp log-likelihood (summed over all spots)
#	Jeremy Silver.
#	21 Jan 2005. Last modified 8 Feb 2005.
{
	beta<-theta[1]
	logsigma<-theta[2]
	logalpha<-theta[3]
	mu <- beta + background
	mu.sf <- foreground - mu - exp(2*logsigma - logalpha)

	dlogdbeta <- -2 * sum(exp(-logalpha) - exp(dnorm(0,mu.sf,exp(logsigma),log = TRUE) - pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE,log.p = TRUE)))
	dlogdlogsigma <- -2 * sum(exp(2*logsigma - 2*logalpha) - (2*exp(2*logsigma - logalpha) + mu.sf)*exp( dnorm(0,mu.sf,exp(logsigma),log = TRUE) - pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE,log.p = TRUE)))
	dlogdlogalpha <- -2 * sum(-1 +(foreground - mu)*exp( -logalpha) - exp(2*logsigma - 2*logalpha) + exp(2*logsigma - logalpha+dnorm(0,mu.sf,exp(logsigma),log = TRUE) - pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE,log.p = TRUE)))

	c(dlogdbeta,dlogdlogsigma,dlogdlogalpha)
}

sumloglik <- function(theta,foreground,background)
#	Minus the norm-exp log-likelihood (summed over all spots)
#	Jeremy Silver.
#	21 Jan 2005. Last modified 8 Feb 2005.
{

	beta<-theta[1]
	logsigma<-theta[2]
	logalpha<-theta[3]
	mu <- beta + background
	mu.sf <- foreground - mu - exp(2*logsigma - logalpha)
	
	-2*sum(-logalpha - (foreground - mu)/exp(logalpha) + 0.5*exp(2*logsigma - 2*logalpha) + log(pnorm(0,mu.sf,exp(logsigma),lower.tail = FALSE)))

}

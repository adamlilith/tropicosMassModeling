#' Estimate number of MCMC samples necessary for convergence
#' 
#' This function is a wrapper for the \code{\link[coda]{raftery.diag}} function in the \pck{coda} package which estimates the number iterations (plus burnin) and thinning interval for MCMC samples necessary for sufficient convergence. This function adds some automation to the process by allowing the user to lower their standards.
#' @param comp A compiled MCMC object produced by \code{\link[nimble]{compileNimble}} of class \code{MCMC_refClass}.
#' @param inits A list of initialization values for the NIMBLE model in \code{mcmc}.
#' @param niter Number of iterations to use for the MCMC used to implement the RL method.
#' @param nburnin Number of burn-in iterations.
#' @param nsamples Number of samples desired in final model (used to determine \code{thin}).
#' @param thin Amount by which to thin MCMC chain used to implement the RL method.
#' @param maxIter Maximum number of iterations to allow in output. The final number of recommended iterations will be \code{min(maxIter, rl)} where \code{rl} is the number recommend by the RL method, rounded up to the nearest 1000. If \code{rl} is > \code{maxIter} a warning will be issued.
#' @param propInsufficient Proportion of coefficients that do not need to be sufficiently sampled when deciding number of iterations. To ensure all are (in theory), this should be 0.
#' @param exact A character list or vector of characters with variable names to ignore when assessing convergence. These names will be matched exactly, so, for example, \code{'p'} will only match a variable named \code{'p'}. If \code{NULL} (default), all monitored variables are used.
#' @param pattern A character list or vector of characters with variable names to ignore when assessing convergence. These names will be matched by pattern, so, for example, \code{'p'} will match \code{'p'}, \code{'p[1]'}, and \code{'psi'}. To ignore series of variables (i.e., with indices), use something like \code{'p[[]'} to remove \code{'p[1]'}, \code{'p[2]'}, \code{'p[3]'}... If \code{NULL} (default), all monitored variables are used.
#' @param retry Logical, if \code{TRUE} (default), and the RF method indicates that the number of samples was insufficient to conduct a reliable RF test, the MCMC is run again with \code{mult} times more samples, repeatedly, until the RF method can be implemented.
#' @param mult Proportion by which to increase number of iterations each attempt.
#' @param verbose If \code{TRUE} display messages. If this is \code{FALSE}, warnings will still be displayed.
#' @param ... arguments to pass to \code{\link[coda]{raftery.diag}}.
#' @return List with these elements:
#' \itemize{
#'		\item{mcmc} Object of class \code{mcmc} with one chain. You can keep this as a chain, if you want to.
#'		\item{sufficient}: \code{TRUE} or \code{FALSE}, indicating if the trial MCMC chain was long enough to implement the RL method.
#'		\item{rlIter}: Number of iterations suggested by the RL method.
#'		\item{recThin}: Thin value if using the \code{rlIter} iterations assuming \code{nsamples} is desired.
#'		\item{rRec}: Number of iterations recommended using the settings in this function.
#' }
#' @examples
#' @export

rafLewis <- function(
	comp,
	inits,
	niter = 2000,
	nburnin = 1000,
	nsamples = 1000,
	thin = 1,
	maxIter = 11000,
	propInsufficient = 0,
	exact = NULL,
	pattern = NULL,
	retry = TRUE,
	mult = 0.2,
	verbose = FALSE,
	...
) {

	if (verbose) omnibus::say('Estimating necessary number of iterations using Raftery-Lewis method...')

	# MCMC
	mcmc <- nimble::runMCMC(
		comp,
		niter = niter,
		nburnin = nburnin,
		thin = thin,
		nchains = 1,
		inits = inits,
		progressBar = verbose,
		samplesAsCodaMCMC = TRUE,
		summary = FALSE,
		WAIC = FALSE,
		...
	)
	flush.console()

	mcmc <- .removeVariablesFromMcmc(mcmc, exact=exact, pattern=pattern)
	mcmcMat <- as.matrix(mcmc)
	
	# RL method
	args <- list(data = mcmcMat, ...)
	dotNames <- omnibus::ellipseNames(...)
	if (any(dotNames == 'q')) args <- c(args, q=dots$q)
	if (any(dotNames == 'r')) args <- c(args, r=dots$r)
	if (any(dotNames == 's')) args <- c(args, s=dots$s)
	if (any(dotNames == 'converge.eps')) args <- c(args, converge.eps=dots$converge.eps)
	
	rafStats <- do.call(coda::raftery.diag, args=args)
	sufficient <- (rafStats$resmatrix[1] != 'Error')

	# repeat if insufficient
	if (!sufficient) {
	
		rlIter <- recIter <- recThin <- NA
		if (!retry) warning('Too few iterations were used for a reliable implementation of the Raftery-Lewis method.')
			
	}
	
	if (!sufficient & retry) {
		
		while (!sufficient & niter <= maxIter) {

			add <- round(mult * niter)
			niter <- add + niter
			# thin <- max(1, floor((niter - nburnin) / nsamples))
			thisThin <- 1
		
			if (verbose) omnibus::say('Increasing iterations to ', niter, ' (before thinning, including burn-in)...')
			mcmc <- addToMcmc(comp = comp, mcmc = mcmc, inits = inits, add = add, thin = thisThin)
			mcmc <- .removeVariablesFromMcmc(mcmc, exact=exact, pattern=pattern)
			args$data <- as.matrix(mcmc)
			rafStats <- do.call(coda::raftery.diag, args=args)
			sufficient <- (rafStats$resmatrix[1] != 'Error')
		
		}
	
	}
	
	if (sufficient) {
		
		rafIters <- rafStats$resmatrix[ , 'N']

		rlIter <- max(rafIters, na.rm=TRUE) + nburnin
		recIter <- quantile(rafIters, probs=1 - propInsufficient, na.rm=TRUE)
		
		recIter <- 1000 * ceiling(recIter / 1000)
		recIter <- max(recIter, nsamples)
		
		if (maxIter < recIter) warning('Number of iterations necessary to achieve desired precision is < maximum desired (maxIter).')
		recIter <- min(maxIter, recIter)
		recIter <- recIter + nburnin
		
		recThin <- max(1, floor((recIter - nburnin) / nsamples))
		
		if (verbose) {
			omnibus::say('Number of iterations + burn-in recommended for adequate sampling ............ ', rlIter)
			omnibus::say('Number of iterations + burn-in recommended based sufficiency ................ ', recIter)
			omnibus::say('Thin interval necessary desired number of samples and degree of precision ... ', recThin)
		}

	}		
		
	out <- list(
		mcmc = mcmc,
		sufficient = sufficient,
		rlIter = rlIter,
		recIter = recIter,
		recThin = recThin
	)
	
	out
	
}

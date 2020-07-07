#' Add iterations to MCMC chains
#'
#' This function is somewhat of a "cheat"... it takes a its main arguments a compiled object of class \code{MCMC_refClass} made by the function  \code{\link[nimble]{compileNimble}} and an object of class \code{\link[coda]{mcmc.list}} and adds iterations to each chain by using the last values in each chain as initialization values. It will not do any burn-in, so adding many iterations could be inefficient.
#' @param comp A compiled MCMC object produced by \code{\link[nimble]{compileNimble}} of class \code{MCMC_refClass}.
#' @param mcmc An object of class \code{\link[coda]{mcmc}} or \code{\link[coda]{mcmc.list}}.
#' @param inits A list object with initialization values. These values are not used; rather, the structure of the initialization values (i.e., names, number of elements) is used to create a new initialization object for each run.
#' @param add Number of additional iterations to add. Note that iterations will be added \emph{after} thinning, so if \code{add} is 1000 and \code{thin} is 10, 100 additional iterations will be added.
#' @param thin Thinning interval.
#' @param verbose If \code{TRUE} display progress. Default is \code{FALSE}.
#' @details Note that adding additional iterations may invalidate and WAIC that was calculated.
#' @return Object of class \code{\link[coda]{mcmc.list}}.
#' @examples
#' @export

addToMcmc <- function(
	comp,
	mcmc,
	inits,
	add = 1000,
	thin = 1,
	verbose = FALSE
) {
	
	existingIters <- coda::niter(mcmc)
	nchain <- coda::nchain(mcmc)
	
	# for each chain
	for (chain in nchains) {

		if (verbose) omnibus::say('Adding to chain ', chain, '...')
	
		# assign new inits values equal to last value in MCMC
		newInits <- inits
		initNames <- names(inits)
		for (i in seq_along(initNames)) {
			thisName <- names(inits)[i]
			initSize <- length(inits[[i]])
			colIndex <- if (initSize == 1) {
				which(colnames(mcmc) == thisName)
			} else {
				match(paste0(thisName, '[', 1:initSize, ']'), colnames(mcmc))
			}
			
			newInits[[thisName]] <- mcmc[existingIters, colIndex]
		}

		mcmcAdd <- nimble::runMCMC(
			mcmc = comp,
			niter = add,
			nburnin = 0,
			thin = 1,
			nchains = 1,
			inits = newInits,
			progressBar = verbose,
			samplesAsCodaMCMC = TRUE,
			summary = FALSE,
			WAIC = FALSE
		)
		flush.console()

		mcmcAdd <- as.matrix(mcmcAdd)
		if (thin > 1) {
			mcmcAdd <- mcmcAdd[seq(1, add, by=thin), ]
		}

		if (class(mcmc) == 'mcmc') {

			mcmc <- rbind(mcmc, mcmcAdd[ , match(colnames(mcmc), colnames(mcmcAdd))])
			mcmc <- coda::as.mcmc(mcmc, thin=thin)
			
		} else if (class(mcmc) == 'mcmc.list') {
		
			mcmc[[chain]] <- rbind(mcmc[[chain]], mcmcAdd[ , match(colnames(mcmc), colnames(mcmcAdd))])
			mcmc <- coda::as.mcmc(mcmc[[chain]], thin=thin)
		
		}

		existingIters <- coda::niter(mcmc)
		
	}
		
	mcmc
		
}


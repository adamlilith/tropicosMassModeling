#' Remove variables from an MCMC chain
#'
#' This function removes variables from an MCMC chain or set of chains.
#' @param mcmc Object of class \code{'mcmc'} or \code{'mcmc.list'}.
#' @param exact A character list or vector of characters with variable names to ignore when assessing convergence. These names will be matched exactly, so, for example, \code{'p'} will only match a variable named \code{'p'}. If \code{NULL} (default), all monitored variables are used.
#' @param pattern A character list or vector of characters with variable names to ignore when assessing convergence. These names will be matched by pattern, so, for example, \code{'p'} will match \code{'p'}, \code{'p[1]'}, and \code{'psi'}. To ignore series of variables (i.e., with indices), use something like \code{'p[[]'} to remove \code{'p[1]'}, \code{'p[2]'}, \code{'p[3]'}... If \code{NULL} (default), all monitored variables are used.
#' @return Object of class \code{'mcmc'} or \code{'mcmc.list'}.
.removeVariablesFromMcmc <- function(
	mcmc,
	exact = NULL,
	pattern = NULL
) {

	if ('mcmc.list' %in% class(mcmc)) {
		for (chain in 1:coda::nchain(mcmc)) {
			mcmc[[chain]] <- .removeVariablesFromMcmc(mcmc[[chain]], exact = exact, pattern = pattern)
		}
	} else {

		# ignore specific variables
		if (!is.null(exact)) {
			for (ign in exact) {
				bads <- which(colnames(mcmc) == ign)
				if (length(bads) > 0) {
					mcmc <- mcmc[ , -bads]
				}
			}
		}

		# ignore patterns
		if (!is.null(pattern)) {
			for (ign in pattern) {
				bads <- which(grepl(colnames(mcmc), pattern=ign))
				if (length(bads) > 0) {
					mcmc <- mcmc[ , -bads]
				}
			}
		}
		
		mcmc <- coda::as.mcmc(mcmc)
		
	}

	mcmc
	
}

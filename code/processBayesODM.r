#' Process output from a trainnBayesODM model
#'
#' Summarizes output for all types of trainBayesODM model functions.
#' @param shape SpatialPolygonsDataFrame.
#' @param mcmc Object of class \code{mcmc.list}.
#' @param effort Character, name of field in \code{shape} with effort.
#' @param detect Character, name of field in \code{shape} with detections.
#' @param hasNeighs Logical vector, one per row in \code{shape}. Indicates if a polygon in \code{shape} has neighbors or not.
#' @param pByState Logical, if \code{FALSE} (default), it is assumed there is one value of p (global) or one value of p per county. If \code{TRUE}, then it is assumed that there is one value of p per state/province, and \code{stateProv} must be specified.
#' @param stateProv Name of field in \code{shape} with state/province names. Ignored if \code{pByState} is \code{FALSE}.
#' @return List.
#' @examples

.processBayesODM <- function(
	shape,
	mcmc,
	effort,
	detect,
	hasNeighs,
	pByState = FALSE,
	stateProv = NULL
) {

	nchains <- length(mcmc$samples)
	
	### remove unneeded stuff
	#########################

		# remove unneeded
		ignores <- c('logit_p[[]', 'z[[]', 'qConstraint')
		for (ignore in ignores) {
			bads <- grepl(colnames(mcmc$samples$chain1), pattern=ignore)
			if (any(bads)) {
				for (chain in 1:nchains) {
					mcmc$samples[[chain]] <- mcmc$samples[[chain]][ , !bads]
				}
			}	
		}
		
		# remove island coefficients for non-islands and CAR coefficients for islands
		if (any(!hasNeighs)) {
		
			# psi_car for islands
			index <- which(!hasNeighs)
			bads <- paste0('psi_car[', index, ']')
			for (chain in 1:nchains) {
				remove <- which(colnames(mcmc$samples[[chain]]) %in% bads)
				mcmc$samples[[chain]] <- mcmc$samples[[chain]][ , -remove]
			}
		
			# psi_island for non-islands
			index <- which(hasNeighs)
			bads <- paste0('psi_island[', index, ']')
			for (chain in 1:nchains) {
				remove <- which(colnames(mcmc$samples[[chain]]) %in% bads)
				mcmc$samples[[chain]] <- mcmc$samples[[chain]][ , -remove]
			}
			
		# no islands
		} else {
		
			# remove psi_island and psi_islandMean for all
			index <- which(hasNeighs)
			bads <- c(paste0('psi_island[', index, ']'), 'psi_islandMean')
			for (chain in 1:nchains) {
				remove <- which(colnames(mcmc$samples[[chain]]) %in% bads)
				mcmc$samples[[chain]] <- mcmc$samples[[chain]][ , -remove]
			}
		
		}

	### calculate summary
	#####################
	
		mcmc$summary$all.chains <- matrix(NA, nrow=ncol(mcmc$samples[[1]]), ncol=6)
		rownames(mcmc$summary$all.chains) <- colnames(mcmc$samples[[1]])
		colnames(mcmc$summary$all.chains) <- c('Mean', 'Median', 'St.Dev', '95%CI_low', '95%CI_upp', 'rhat')
		
		samps <- mcmc$samples[[1]]
		if (nchains > 1) {
			for (chain in 2:nchains) {
				samps <- rbind(samps, mcmc$samples[[chain]])
			}
		}
		
		mcmc$summary$all.chains[ , 'Mean'] <- colMeans(samps, na.rm=TRUE)
		mcmc$summary$all.chains[ , 'Median'] <- apply(samps, 2, median, na.rm=TRUE)
		mcmc$summary$all.chains[ , 'St.Dev'] <- apply(samps, 2, sd, na.rm=TRUE)
		mcmc$summary$all.chains[ , '95%CI_low'] <- apply(samps, 2, quantile, probs=0.025, na.rm=TRUE)
		mcmc$summary$all.chains[ , '95%CI_upp'] <- apply(samps, 2, quantile, probs=0.975, na.rm=TRUE)
		# mcmc$summary$all.chains[ , 'rhat'] <- tryCatch(
			# expr={wiqid::simpleRhat(mcmc$samples, n.chains=nchains)},
			# error = function(x) NA,
			# finally={rep(NA, nrow(mcmc$summary$all.chains))}
		# )

		rm(samps); gc()

	### update shapefile
	####################

		# psi
		psiIndex <- which(grepl(rownames(mcmc$summary$all.chains), pattern='psi[[]'))
		psi <- mcmc$summary$all.chains[psiIndex, 'Mean']
		psiUncer <- mcmc$summary$all.chains[psiIndex, '95%CI_upp'] - mcmc$summary$all.chains[psiIndex, '95%CI_low']
		
		# p
		pIndex <- if (any(grepl(rownames(mcmc$summary$all.chains), pattern='p[[]'))) {
			which(grepl(rownames(mcmc$summary$all.chains), pattern='p[[]'))
		} else {
			which(rownames(mcmc$summary$all.chains) == 'p')
		}
		p <- mcmc$summary$all.chains[pIndex, 'Mean']
		pUncer <- mcmc$summary$all.chains[pIndex, '95%CI_upp'] - mcmc$summary$all.chains[pIndex, '95%CI_low']
		
		# assume p values pertain to states, not counties
		if (pByState) {
		
			state <- as.numeric(as.factor(shape@data[ , stateProv]))
			p <- p[state]
			pUncer <- pUncer[state]

		}
		
		shape@data$psi <- psi
		shape@data$psiUncer <- psiUncer
		shape@data$p <- p
		shape@data$pUncer <- pUncer
		
		names(shape@data)[(ncol(shape@data) - 3):ncol(shape@data)] <- c('psi', 'psi90CI', 'p', 'p90CI')
		
	### done
	########
	
		list(mcmc=mcmc, shape=shape)

}

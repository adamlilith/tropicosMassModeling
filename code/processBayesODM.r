#' Process output from a trainnBayesODM model
#'
#' Summarizes output for all types of trainBayesODM model functions.
#' @param shape SpatialPolygonsDataFrame.
#' @param mcmc Object of class \code{mcmc.list}.
#' @param effort Character, name of field in \code{shape} with effort.
#' @param detect Character, name of field in \code{shape} with detections.
#' @param neighIndex Index of units with neighbors (from the spatial object). Set as \code{NULL} if there are none.
#' @param islandndex Index of units that are islands (from the spatial object). Set as \code{NULL} if there are none.
#' @param pByState Logical, if \code{FALSE} (default), it is assumed there is one value of p (global) or one value of p per county. If \code{TRUE}, then it is assumed that there is one value of p per state/province, and \code{stateProv} must be specified.
#' @param stateProv Name of field in \code{shape} with state/province names. Ignored if \code{pByState} is \code{FALSE}.
#' @return List.
#' @examples

.processBayesODM <- function(
	shape,
	mcmc,
	effort,
	detect,
	neighIndex = NULL,
	islandIndex = NULL,
	pByState = FALSE,
	stateProv = NULL
) {

	nchains <- coda::nchain(mcmc)
	mcmc <- list(samples=mcmc)
	
	### calculate summary
	#####################
	
		mcmc$summary$all.chains <- matrix(NA, nrow=ncol(mcmc$samples[[1]]), ncol=5)
		rownames(mcmc$summary$all.chains) <- colnames(mcmc$samples[[1]])
		colnames(mcmc$summary$all.chains) <- c('Mean', 'Median', 'St.Dev', '95%CI_low', '95%CI_upp')
		
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

		rm(samps); gc()

	### update shapefile
	####################

		# psi
		if (!is.null(neighIndex[1])) {
		
			psiNeighIndex <- which(grepl(rownames(mcmc$summary$all.chains), pattern='psi_neigh[[]'))
			psiNeigh <- mcmc$summary$all.chains[psiNeighIndex, 'Mean']
			psiUncerNeigh <- mcmc$summary$all.chains[psiNeighIndex, '95%CI_upp'] - mcmc$summary$all.chains[psiNeighIndex, '95%CI_low']
		
		}
		
		if (!is.null(islandIndex[1])) {
			
			psiIslandIndex <- which(grepl(rownames(mcmc$summary$all.chains), pattern='psi_island[[]'))
			psiIsland <- mcmc$summary$all.chains[psiIslandIndex, 'Mean']
			psiUncerIsland <- mcmc$summary$all.chains[psiIslandIndex, '95%CI_upp'] - mcmc$summary$all.chains[psiIslandIndex, '95%CI_low']
		
		}
			
		# p
		pIndex <- if (any(grepl(rownames(mcmc$summary$all.chains), pattern='p[[]'))) {
			which(grepl(rownames(mcmc$summary$all.chains), pattern='p[[]'))
		} else {
			which(rownames(mcmc$summary$all.chains) == 'p')
		}
		p <- mcmc$summary$all.chains[pIndex, 'Mean']
		pUncer <- mcmc$summary$all.chains[pIndex, '95%CI_upp'] - mcmc$summary$all.chains[pIndex, '95%CI_low']
		
		# if p values pertain to states, not counties
		if (pByState) {
		
			state <- as.numeric(as.factor(shape@data[ , stateProv]))
			p <- p[state]
			pUncer <- pUncer[state]

		}
		
		# assign
		shape@data$pUncer <- shape@data$p <- shape@data$psiUncer <- shape@data$psi <- NA
		
		shape@data$psiUncer <- shape@data$psi <- NA
		if (!is.null(neighIndex[1])) {
			shape@data$psi[neighIndex] <- psiNeigh
			shape@data$psiUncer[neighIndex] <- psiUncerNeigh
		}
		
		if (!is.null(islandIndex[1])) {
			shape@data$psi[islandIndex] <- psiIsland
			shape@data$psiUncer[islandIndex] <- psiUncerIsland
		}
		
		shape@data$p <- p
		shape@data$pUncer <- pUncer
		
		names(shape@data)[(ncol(shape@data) - 3):ncol(shape@data)] <- c('psi', 'psi90CI', 'p', 'p90CI')
		
	### done
	########
	
		list(mcmc=mcmc, shape=shape)

}

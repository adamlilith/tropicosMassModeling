#' Calculate log likelihood for each observation of an enmSdmBayes model
#'
#' This function calculates the log likelihood of each observation (detections within a county).
#' @param mcmc Object of class \code{mcmc.list}.
#' @param data Data list object (like that supplied to \code{nimbleMCMC} or \code{runMCMC}).
#' @param constants Constants list object (like that supplied to \code{nimbleMCMC} or \code{runMCMC}).
#' @param hasNeighs Logical, if \code{TRUE} then at least two counties are spatial neighbors.
#' @param hasIslands Logical, if \code{TRUE} then at least one county has no spatial neighbors.
#' @expor

likeFromMcmc <- function(
	mcmc,
	data,
	constants,
	hasNeighs = TRUE,
	hasIslands = TRUE
) {

	# combine chains
	nchains <- coda::nchain(mcmc)
	niter <- coda::niter(mcmc) * nchains
	combined <- mcmc[[1]]
	if (nchains > 1) {
		for (chain in 2:nchains) {
			combined <- rbind(combined, mcmc[[chain]])
		}
	}
	
	# counties with neighborhoods (not islands)
	if (hasNeighs) {
	
		neigh_ll <- matrix(NA, nrow=niter, ncol=length(data$y_neigh))

		# detection and occupancy
		for (iter in 1:niter) {

			# CAR component
			psi_tau <- combined[iter, 'psi_tau', drop=TRUE]
			psi_tau_ll <- dgamma(psi_tau, 0.001, 0.001, log=TRUE)
			
			psi_car <- combined[iter, grepl(colnames(combined), pattern='psi_car[[]')]
			psi_car_ll <- dcar_normal(psi_car, adj=data$carAdjCounty[1:constants$lengthAdjCounties], weights=data$carWeightCounty[1:constants$lengthAdjCounties], num=data$carNumCounty[1:constants$numNeighs], tau=psi_tau, c=3, zero_mean=0, log=TRUE)

			# detection
			z_neigh <- combined[iter, grepl(colnames(mcmc), pattern='z_neigh[[]'), drop=TRUE]
			p <- combined[iter, grepl(colnames(mcmc), pattern='p[[]'), drop=TRUE]
			q <- combined[iter, 'q', drop=TRUE]

			p_ll <- dbeta(p, 1, 1)
			q_ll <- dbeta(q, 1, 2)
			
			p_star_neigh <- z_neigh * p[constants$state_neigh] + (1 - z_neigh[g]) * q
			
			neigh_ll[iter, ] <-
				dbin(data$y_neigh, p_star_neigh, data$N_neigh, log = TRUE) + 			# detection
				p_ll[constants$state_neigh] +									# state effect
				q_ll +															# false detection
				psi_car															# occupancy ~ CAR

		} # next iteration
		
	} # has neighbors

	# island counties
	if (hasIslands) {
	
		island_ll <- matrix(NA, nrow=niter, ncol=length(data$y_island))

		# detection and occupancy
		for (iter in 1:niter) {

			# detection
			z_island <- combined[iter, grepl(colnames(mcmc), pattern='z_island[[]'), drop=TRUE]
			p <- combined[iter, grepl(colnames(mcmc), pattern='p[[]'), drop=TRUE]
			q <- combined[iter, 'q', drop=TRUE]

			p_ll <- dbeta(p, 1, 1)
			q_ll <- dbeta(q, 1, 2)
			
			p_star_island <- z_island[h] * p[constants$state_island[h]] + (1 - z_island[h]) * q
			
			psi_islandMean <- combined[iter, 'psi_islandMean', drop=TRUE]
			psi_islandSd <- combined[iter, 'psi_islandSd', drop=TRUE]
			
			psi_islandMean_ll <- dnorm(psi_islandMean, mean=0, sd=1000)
			psi_islandMean_ll <- dgamma(psi_islandSd, 0.001, 0.001)
			
			island_ll[iter, ] <-
				dbin(data$y_island, p_star_island, data$N_island, log = TRUE) +	# detection
				p_ll[constants$state_island] +									# state effect
				q_ll +															# false detection
				psi_islandMean_ll +												# psi for islands
				psi_islandSd_ll													# psi for islands

		} # next iteration
		
	} # if has islands

	# return
	if (hasNeighs & !hasIslands) {
		neigh_ll
	} else if (!hasNeighs & hasIslands) {
		island_ll
	} else {
		cbind(neigh_ll, island_ll)
	}
	
}
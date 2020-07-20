#' Create a bias-corrected distribution model
#' 
#' This function implements an occupancy-detection model assuming that:
#' \itemize{
#' 		\item{Detection varies by state but not between counties in a state.
#' 		\item{False detections can occur when the species does not occupy a county but does not occur when it does.
#'		\item{Occupancy is a function of the occupancy of neighboring counties.
#' }
#' @param shape SpatialPolygonsDataFrame
#' @param effort Character, name of field in \code{shape} indicating collection effort.
#' @param detect Character, name of field in \code{shape} indicating number of detections of focal species.
#' @param stateProv Character, name of field with state/province.
#' @param county Character, name of field with county.
#' @param niter Positive integer, number of MCMC iterations (including burn-in). The default is 2000, but this is often too low for most cases (i.e., try numbers in the range 10000, 100000, etc.).
#' @param nburnin Positive integer, number of burn-in samples (less than \code{niter}). The default is 1000, but this is often too low for most cases (i.e., try numbers in the range 10000, 100000, etc.).
#' @param nchains Positive integer, number of MCMC chains (default is 4).
#' @param thin Positive integer, number of MCMC samples by which to thin the results (default is 1; i.e., no thinning). To reduce memory requirements, you can use \code{thin} while increasing \code{niter} and/or \code{nburnin}.
#' @param raf Logical, if \code{TRUE} (default), dynamically change the value for niter using the Raftery-Lewis method for assessing number of iterations needed for convergence. The number used will the the number suggested rounded up to the nearest 1000. The total number of iterations (including burn-in) will be this number plus \code{nburnin}.
#' @param minIter Minimum number of iterations to require (not including burn-in).
#' @param maxIter Maximum number of iterations to allow in output. The final number of recommended iterations will be \code{min(maxIter, rl)} where \code{rl} is the number recommend by the RL method, rounded up to the nearest 1000. If \code{rl} is > \code{maxIter} a warning will be issued.
#' @param nsamples Number of samples desired (after burn-in and thinning). Only used if \code{raf} is \code{TRUE}. Used to reset the value of \code{thin} after estimating the number of iterations using the Raftery-Lewis method.
#' @param na.rm Logical, if \code{TRUE} (default), then remove rows in \code{shape} that have \code{NA} in any input field. If this is \code{FALSE} and \code{NA}s occur in any fields then an error may occur.
#' @param verbose Logical, if \code{TRUE} (default) display progress.
#' @param ... Arguments to pass to \code{\link[nimble]{configureMCMC}}, \code{\link[nimble]{runMCMC}}, and \code{rafLewis}.
#' @return A list object with three elements, one with the MCMC chains, a second with the object \code{shape} with model output appended to the data portion of the object, and a third with model information. The new columns in \code{shape} represent:
#' \itemize{
#' 		\item{\code{psi}} Probability of occurrence
#' 		\item{\code{psi95CI}} 95% credibility interval for probability of occurrence
#' 		\item{\code{p}} Probability of detecting the focal species assuming it is present
#' 		\item{\code{p95CI}} 95% credibility interval for probability of detection
#' }
#' @examples
#' @export

trainBayesODM_pStateMean_psiCar <- function(
	shape,
	effort,
	detect,
	stateProv,
	county,
	niter = 2000,
	nburnin = 1000,
	nchains = 4,
	thin = 1,
	raf = TRUE,
	minIter = 2000,
	maxIter = 11000,
	nsamples = 1000,
	na.rm = TRUE,
	verbose = TRUE,
	...
) {

	### prepare data
	################

		# remove missing
		if (na.rm) {
			ok <- complete.cases(shape@data[ , c(effort, detect)])
			shape <- shape[ok, ]
		}
			
		### CAR setup for counties
		
		neighs <- spdep::poly2nb(shape, queen=TRUE)
		neighIndex <- which(!sapply(neighs, function(shape) { length(shape == 1) && (shape == 0) }))
		islandIndex <- which(sapply(neighs, function(shape) { length(shape == 1) && (shape == 0) }))

		hasNeighs <- length(neighIndex) > 0
		hasIslands <- length(islandIndex) > 0
		
		if (hasIslands) {
			shapeNoIslands <- shape[neighIndex, ]
			neighs <- spdep::poly2nb(shapeNoIslands, queen=TRUE)
		}

		numNeighs <- length(neighIndex)
		numIslands <- length(islandIndex)
		
		numCounties <- nrow(shape)
		numStates <- length(unique(shape@data[ , stateProv]))
		
		### input
		data <- list()
		constants <- list(
			numStates = numStates
		)
		
		inits <- list(
			q = 0.01,
			qConstraint = rep(1, numStates),
			p = runif(numStates, 0.1, 0.3)
		)

		# some counties have neighbors
		if (hasNeighs) {
			
			# remove islands
			neighsList <- spdep::nb2WB(neighs)
			carAdjCounty <- neighsList$adj
			carWeightCounty <- neighsList$weights
			carNumCounty <- neighsList$num

			data <- c(
				data,
				list(
					y_neigh = shape@data[neighIndex, detect],
					N_neigh = shape@data[neighIndex, effort],
					carAdjCounty = carAdjCounty,
					carWeightCounty = carWeightCounty,
					carNumCounty = carNumCounty
				)
			)

			constants <- c(
				constants,
				list(
					numNeighs = numNeighs,
					state_neigh = as.numeric(as.factor(shape@data[neighIndex, stateProv])),

					lengthAdjCounties = length(carAdjCounty)
				)
			)

			inits <- c(
				inits,
				list(
					p_star_neigh = runif(numNeighs, 0.1, 0.5),
				
					z_neigh = rep(1, numNeighs),
					psi_neigh = runif(numNeighs, 0.4, 0.6),
					
					psi_tau = 0.1,
					psi_car = rnorm(numNeighs)
				)
			)
			
		}
		
		# some counties are islands
		if (hasIslands) {
		
			data <- c(
				data,
				list(
					y_island = shape@data[islandIndex, detect],
					N_island = shape@data[islandIndex, effort]
				)
			)

			numStates_island <- length(shape@data[islandIndex, stateProv])
			
			constants <- c(
				constants,
				list(
					numIslands = numIslands,
					state_island = as.numeric(as.factor(shape@data[islandIndex, stateProv]))
				)
			)
			
			psi_island <- runif(numIslands, 0.4, 0.6)
			
			inits <- c(
				inits,
				list(
					p_star_island = runif(numIslands, 0.1, 0.5),
					z_island = rep(1, numIslands),
					psi_island = psi_island,
					logit_psi_island = logit(psi_island),
					psi_islandMean = 0,
					psi_islandTau = 1
				)
			)
			
		}
		
		### monitors
		monitors <- names(inits)
		# monitors <- monitors[!(monitors %in% c('z_neigh', 'z_neigh'))]
		monitors <- monitors[!(monitors %in% c('logit_psi_island'))]

	### model setup
	###############
	
		code <- if (hasNeighs & hasIslands) {
			.trainBayesODM_pStateMean_psiCar_neighsIslands
		} else if (hasNeighs & !hasIslands) {
			.trainBayesODM_pStateMean_psiCar_neighsOnly
		} else if (!hasNeighs & hasIslands) {
			.trainBayesODM_pStateMean_psiCar_islandsOnly
		}
		
		### construct and compile model
		###############################
		
			model <- nimble::nimbleModel(code=code, constants=constants, data=data, inits=inits, check=TRUE)
			flush.console()
		
			conf <- nimble::configureMCMC(model, monitors = monitors, print = verbose, enableWAIC = TRUE)
			
			### modify samplers
			node <- 'q'
			conf$removeSamplers(node)
			conf$addSampler(target=node, type='slice')
			
			if (hasNeighs) {
				node <- 'psi_tau'
				conf$removeSamplers(node)
				conf$addSampler(target=node, type='slice')
			}
			
			if (hasIslands) {
				node <- 'psi_islandTau'
				conf$removeSamplers(node)
				conf$addSampler(target=node, type='slice')

				# for (i in 1:numIslands) {
					# node <- paste0('logit_psi_island[', i, ']')
					# conf$removeSamplers(node)
					# conf$addSampler(target=node, type='slice')
				# }
			}
			
			for (i in 1:numStates) {

				node <- paste0('p[', i, ']')
				conf$removeSamplers(node)
				conf$addSampler(target=node, type='slice')
				
			}
			
			confBuild <- nimble::buildMCMC(conf)
			compiled <- nimble::compileNimble(model, confBuild)

	### variables to ignore when assessing convergence
	##################################################

		exact <- NULL
		# pattern <- c('qConstraint[[]', 'p_star_neigh[[]', 'z_neigh[[]', 'p_star_island[[]', 'z_island[[]', 'logit_psi_island[[]')
		pattern <- c('qConstraint[[]', 'p_star_neigh[[]', 'p_star_island[[]', 'logit_psi_island[[]')
	
	### use Raftery-Lewis statistic to see how many iterations are needed
	#####################################################################

		if (raf) {

			rl <- rafLewis(
				comp=compiled$confBuild,
				inits = inits,
				niter = niter,
				nburnin = nburnin,
				nsamples = nsamples,
				thin = thin,
				minIter = minIter,
				maxIter = maxIter,
				retry = TRUE,
				exact = exact,
				pattern = pattern,
				propInsufficient = minUnconverged,
				verbose = verbose
			)

			rl$mcmc <- .removeVariablesFromMcmc(rl$mcmc, exact=exact, pattern=pattern)

			rand <- round(1E6 * runif(1))
			tempFile <- paste0('./temp/temp', rand, '.rda')
			save(rl, file=tempFile)
			rm(rl)
			gc()
			load(tempFile)
			file.remove(tempFile)
			
		}
		
		if (rl$sufficient) {
		
			niter <- rl$recIter
			thin <- rl$recThin
			
			if (coda::niter(rl$mcmc) < nsamples) {
				if (verbose) omnibus::say('Adding iterations to Raftery-Lewis MCMC chain...')
				add <- (niter - burnin - coda::niter(rl$mcmc)) * thin
				rl$mcmc <- addToMcmc(comp = comp, mcmc = rl$mcmc, inits = inits, add = add, thin = thin)
			}
			
		}
			
	### MCMC
	########
		
		if (nchains > 1 | !raf) {

			startChain <- if (raf) { 2 } else { 1 }
			thisNumChains <- if (raf) { nchains - 1 } else { nchains }
		
			# note: thinning afterward
			mcmc <- runMCMC(compiled$confBuild, niter = niter, nburnin = nburnin, thin = thin, nchains = thisNumChains, inits = inits, progressBar = verbose, samplesAsCodaMCMC = TRUE, summary = FALSE, WAIC = FALSE)
			flush.console()
			
			mcmc <- .removeVariablesFromMcmc(mcmc, exact=exact, pattern=pattern)
							
		}
		
		# combine all chains
		nrowRafMcmc <- coda::niter(rl$mcmc)
		nrowNewMcmc <- coda::niter(mcmc)
		
		if (nrowRafMcmc > nrowNewMcmc) {
			rl$mcmc <- rl$mcmc[1:nrowNewMcmc, ]
			rl$mcmc <- as.mcmc(rl$mcmc)
		}
		
		mcmc[[1 + coda::nchain(mcmc)]] <- rl$mcmc
		names(mcmc)[coda::nchain(mcmc)] <- paste0('chain', coda::nchain(mcmc))

	### WAIC and LOO
	################
	
		# combine chains
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
				z_neigh <- combined[iter, grepl(colnames(combined), pattern='z_neigh[[]'), drop=TRUE]
				p <- combined[iter, grepl(colnames(combined), pattern='p[[]'), drop=TRUE]
				q <- combined[iter, 'q', drop=TRUE]

				p_ll <- dbeta(p, 1, 1, log=TRUE)
				q_ll <- dbeta(q, 1, 2, log=TRUE)
				
				p_star_neigh <- z_neigh * p[constants$state_neigh] + (1 - z_neigh) * q
				
				neigh_ll[iter, ] <-
					dbinom(data$y_neigh, size=data$N_neigh, prob=p_star_neigh, log = TRUE) +	# detection
					p_ll[constants$state_neigh] +									# state effect
					q_ll +															# false detection
					psi_car_ll														# occupancy ~ CAR

			} # next iteration
			
		} # has neighbors

		# island counties
		if (hasIslands) {
		
			island_ll <- matrix(NA, nrow=niter, ncol=length(data$y_island))

			# detection and occupancy
			for (iter in 1:niter) {

				# detection
				z_island <- combined[iter, grepl(colnames(combined), pattern='z_island[[]'), drop=TRUE]
				p <- combined[iter, grepl(colnames(combined), pattern='p[[]'), drop=TRUE]
				q <- combined[iter, 'q', drop=TRUE]

				p_ll <- dbeta(p, 1, 1, log=TRUE)
				q_ll <- dbeta(q, 1, 2, log=TRUE)
				
				p_star_island <- z_island * p[constants$state_island] + (1 - z_island) * q
				
				psi_islandMean <- combined[iter, 'psi_islandMean', drop=TRUE]
				psi_islandTau <- combined[iter, 'psi_islandTau', drop=TRUE]
				
				psi_islandMean_ll <- dnorm(psi_islandMean, mean=0, sd=1000, log=TRUE)
				psi_islandSd_ll <- dgamma(psi_islandTau, 1, 0.001, log=TRUE)
				
				island_ll[iter, ] <-
					dbinom(data$y_island, size=data$N_island, prob=p_star_island, log = TRUE) +	# detection
					p_ll[constants$state_island] +									# state effect
					q_ll +															# false detection
					psi_islandMean_ll +												# psi for islands
					psi_islandSd_ll													# psi for islands

			} # next iteration
			
		} # if has islands

		# LL of neighbors and islands, one row per iteration, one column per county
		ll <- if (hasNeighs & !hasIslands) {
			neigh_ll
		} else if (!hasNeighs & hasIslands) {
			island_ll
		} else {
			cbind(neigh_ll, island_ll)
		}
	
		waic <- loo::waic(ll)
		loo <- loo::loo(ll, cores=4)
	
		
	### process model output
	########################

		mcmc <- .removeVariablesFromMcmc(mcmc, exact=exact, pattern=c('z_neigh[[]', 'z_island[[]'))

		mcmcModel <- .processBayesODM(shape=shape, mcmc=mcmc, effort=effort, detect=detect, neighIndex=neighIndex, islandIndex=islandIndex, pByState=TRUE, stateProv=stateProv)

		meta <- list(
			niter=niter,
			nburnin=nburnin,
			nchains=nchains,
			thin=thin,
			effort=effort,
			detect=detect,
			stateProv=stateProv,
			county=county,
			hasIslands=hasIslands,
			na.rm=na.rm,
			...
		)
		
		mcmcModel <- c(mcmcModel, list(code=code), meta=list(meta), waic=list(waic), loo=list(loo))
		mcmcModel

}

.trainBayesODM_pStateMean_psiCar_neighsIslands <- nimble::nimbleCode({
	
	### likelihood for counties with neighborhoods
	for (g in 1:numNeighs) {

		# detection
		y_neigh[g] ~ dbin(p_star_neigh[g], N_neigh[g])
		p_star_neigh[g] <- z_neigh[g] * p[state_neigh[g]] + (1 - z_neigh[g]) * q
		
		# occupancy
		z_neigh[g] ~ dbern(psi_neigh[g])
		logit(psi_neigh[g]) <- psi_car[g]
	}
		
	# county occupancy CAR
	psi_tau ~ dgamma(0.001, 0.001)
	psi_car[1:numNeighs] ~ dcar_normal(adj=carAdjCounty[1:lengthAdjCounties], weights=carWeightCounty[1:lengthAdjCounties], num=carNumCounty[1:numNeighs], tau=psi_tau, c=3, zero_mean=0)
	
	### likelihood for island counties
	for (h in 1:numIslands) {

		# detection
		y_island[h] ~ dbin(p_star_island[h], N_island[h])
		p_star_island[h] <- z_island[h] * p[state_island[h]] + (1 - z_island[h]) * q
		
		# occupancy
		z_island[h] ~ dbern(psi_island[h])
		logit(psi_island[h]) ~ dnorm(psi_islandMean, sd=psi_islandTau)

	}

	psi_islandMean ~ dnorm(0, tau=0.001)
	psi_islandTau ~ dgamma(1, 0.001)
	
	### general priors
	q ~ dbeta(1, 2)
	
	# detection
	for (i in 1:numStates) {
		p[i] ~ dbeta(1, 1)
		qConstraint[i] ~ dconstraint(p[i] > q)
	}
	
})

.trainBayesODM_pStateMean_psiCar_neighsOnly <- nimble::nimbleCode({
	
	### likelihood for counties with neighborhoods
	for (g in 1:numNeighs) {

		# detection
		y_neigh[g] ~ dbin(p_star_neigh[g], N_neigh[g])
		p_star_neigh[g] <- z_neigh[g] * p[state_neigh[g]] + (1 - z_neigh[g]) * q
		
		# occupancy
		z_neigh[g] ~ dbern(psi_neigh[g])
		logit(psi_neigh[g]) <- psi_car[g]
	}
		
	# county occupancy CAR
	psi_tau ~ dgamma(0.001, 0.001)
	psi_car[1:numNeighs] ~ dcar_normal(adj=carAdjCounty[1:lengthAdjCounties], weights=carWeightCounty[1:lengthAdjCounties], num=carNumCounty[1:numNeighs], tau=psi_tau, c=3, zero_mean=0)
	
	### general priors
	q ~ dbeta(1, 2)
	
	# detection
	for (i in 1:numStates) {
		p[i] ~ dbeta(1, 1)
		qConstraint[i] ~ dconstraint(p[i] > q)
	}
	
})

.trainBayesODM_pStateMean_psiCar_islandsOnly <- nimble::nimbleCode({
	
	### likelihood for island counties
	for (h in 1:numIslands) {

		# detection
		y_island[h] ~ dbin(p_star_island[h], N_island[h])
		p_star_island[h] <- z_island[h] * p[state_island[h]] + (1 - z_island[h]) * q
		
		# occupancy
		z_island[h] ~ dbern(psi_island[h])
		logit(psi_island[h]) ~ dnorm(psi_islandMean, sd=psi_islandTau)

	}

	psi_islandMean ~ dnorm(0, tau=0.001)
	psi_islandTau ~ dgamma(1, 0.001)
	
	### general priors
	q ~ dbeta(1, 2)
	
	# detection
	for (i in 1:numStates) {
		p[i] ~ dbeta(1, 1)
		qConstraint[i] ~ dconstraint(p[i] > q)
	}
	
})

#' Create a bias-corrected distribution model
#' 
#' This function implements an occupancy-detection model assuming that:
#' \itemize{
#' 		\item{Detection is constant across counties.
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
#' @param na.rm Logical, if \code{TRUE} (default), then remove rows in \code{shape} that have \code{NA} in any input field. If this is \code{FALSE} and \code{NA}s occur in any fields then an error may occur.
#' @param verbose Logical, if \code{TRUE} (default) display progress.
#' @param ... Arguments to pass to \code{\link[nimble]{configureMCMC}} and \code{\link[nimble]{runMCMC}}.
#' @return A list object with three elements, one with the MCMC chains, a second with the object \code{shape} with model output appended to the data portion of the object, and a third with model information. The new columns in \code{shape} represent:
#' \itemize{
#' 		\item{\code{psi}} Probability of occurrence
#' 		\item{\code{psi95CI}} 95% credibility interval for probability of occurrence
#' 		\item{\code{p}} Probability of detecting the focal species assuming it is present
#' 		\item{\code{p95CI}} 95% credibility interval for probability of detection
#' }
#' @examples
#' @export

trainBayesODM_pConstant_psiCar <- function(
	shape,
	effort,
	detect,
	stateProv,
	county,
	niter = 2000,
	nburnin = 1000,
	nchains = 4,
	thin = 1,
	na.rm=TRUE,
	verbose=TRUE,
	...
) {

	### prepare data
	################

		# remove missings
		if (na.rm) {
			ok <- complete.cases(shape@data[ , c(effort, detect, stateProv, county)])
			shape <- shape[ok, ]
		}
			
		### CAR setup for counties
		
		neighs <- spdep::poly2nb(shape, queen=TRUE)
		
		hasNeighs <- !sapply(neighs, function(shape) { length(shape == 1) && (shape == 0) })
		
		neighsList <- spdep::nb2WB(neighs)
		carAdjCounty <- neighsList$adj
		carWeightCounty <- neighsList$weights
		carNumCounty <- neighsList$num

		### data

		data <- list(
			y = shape@data[ , detect],
			N = shape@data[ , effort],
			carAdjCounty = carAdjCounty,
			carWeightCounty = carWeightCounty,
			carNumCounty = carNumCounty
		)
			
		### constants
		constants <- list(
			yConst = data$y,
			hasNeighs = as.numeric(hasNeighs),
			numCounties = nrow(shape),
			lengthAdjCounties = length(carAdjCounty)
		)
		
		### initializations
		# z <- as.numeric(data$y > 0)
		psi <- runif(constants$numCounties, 0.4, 0.6)
		
		inits <- list(
		
			z = rep(1, constants$numCounties),
			psi = psi,
			
			q = 0.01,
			qConstraint = 1,
			
			# p = runif(constants$numCounties, 0.3, 0.5),
			p_star = 0.7,
			
			psi_island = rnorm(constants$numCounties),
			psi_islandMean = -0.1,
			
			psi_tau = 0.1,
			psi_car = rnorm(constants$numCounties)
			
		)

		### monitors
		monitors <- names(inits)
		monitors <- monitors[!(monitors %in% c('logit_p'))]
		monitors <- monitors[!(monitors %in% c('z'))]

	### model setup
	###############
	
		code <- nimble::nimbleCode({
			
			### likelihood
			for (i in 1:numCounties) {

				# detection
				y[i] ~ dbin(p[i], N[i])
				p[i] <- z[i] * p_star + (1 - z[i]) * q

				# occupancy
				z[i] ~ dbern(psi[i])
				logit(psi[i]) <- hasNeighs[i] * psi_car[i] + (1 - hasNeighs[i]) * psi_island[i]
				psi_island[i] ~ dnorm(psi_islandMean, tau=0.001)
			
			}

			### priors
			q ~ dbeta(1, 10)
			qConstraint ~ dconstraint(q < p_star)
			
			# detection
			p_star ~ dbeta(1, 1)
			
			# occupancy
			psi_islandMean ~ dnorm(0, tau=0.001)
			
			# county occupancy CAR
			psi_tau ~ dgamma(0.001, 0.001)
			psi_car[1:numCounties] ~ dcar_normal(adj=carAdjCounty[1:lengthAdjCounties], weights=carWeightCounty[1:lengthAdjCounties], num=carNumCounty[1:numCounties], tau=psi_tau, c=3, zero_mean=0)
			
		})
		
		
		### construct and compile model
		model <- nimble::nimbleModel(code=code, constants=constants, data=data, inits=inits, check=TRUE)
		flush.console()
		
		# conf <- nimble::configureMCMC(model, monitors = monitors, thin = thin, print = verbose, ...)
		conf <- nimble::configureMCMC(model, monitors = monitors, thin = thin, print = verbose, enableWAIC = TRUE)
		flush.console()
		
		### modify samplers
		node <- 'q'
		conf$removeSamplers(node)
		conf$addSampler(target=node, type='slice')

		node <- 'p_star'
		conf$removeSamplers(node)
		conf$addSampler(target=node, type='slice')

		node <- 'psi_tau'
		conf$removeSamplers(node)
		conf$addSampler(target=node, type='slice')
		
		for (i in 1:constants$numCounties) {
			node <- paste0('psi_island[', i, ']')
			conf$removeSamplers(node)
			conf$addSampler(target=node, type='slice')
		}
		
		confBuild <- nimble::buildMCMC(conf)
		flush.console()
		compiled <- nimble::compileNimble(model, confBuild)
		flush.console()

	### MCMC
	########
		
		# mcmc <- runMCMC(compiled$confBuild, niter = niter, nburnin = nburnin, nchains = nchains, inits = inits, progressBar = verbose, samplesAsCodaMCMC = TRUE, summary = FALSE, ...)
		mcmc <- runMCMC(compiled$confBuild, niter = niter, nburnin = nburnin, nchains = nchains, inits = inits, progressBar = verbose, samplesAsCodaMCMC = TRUE, summary = FALSE, WAIC = TRUE)
		flush.console()

	### process model output
	########################

		mcmcModel <- .processBayesODM(shape=shape, mcmc=mcmc, effort=effort, detect=detect, hasNeighs=hasNeighs)

		meta <- list(
			niter=niter,
			nburnin=nburnin,
			nchains=nchains,
			thin=thin,
			effort=effort,
			detect=detect,
			stateProv=stateProv,
			county=county,
			hasIslands=any(!hasNeighs),
			na.rm=na.rm
		)
		
		mcmcModel <- c(mcmcModel, list(code=code), meta=meta)
		mcmcModel

}

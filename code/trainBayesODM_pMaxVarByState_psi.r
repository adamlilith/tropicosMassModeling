#' Create a bias-corrected distribution model
#' 
#' This function implements an occupancy-detection model assuming that:
#' \itemize{
#' 		\item{Maximum detectability varies by state.
#' 		\item{False detections can occur when the species does not occupy a county but does not occur when it does.
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

trainBayesODM_pMaxVarByState_psi <- function(
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
			
		### data

		data <- list(
			y = shape@data[ , detect],
			N = shape@data[ , effort]
		)
			
		### constants
		constants <- list(
			yConst = data$y,
			numCounties = nrow(shape),
			numStates = length(unique(shape@data[ , stateProv])),
			state = as.numeric(as.factor(shape@data[ , stateProv]))
		)
		
		### initializations
		psi <- runif(constants$numCounties, 0.4, 0.6)
		p <- runif(constants$numCounties, 0.1, 0.3)
		
		inits <- list(
		
			z = rep(1, constants$numCounties),
			psi = psi,
			p = p,
			
			q = 0.01,
			qConstraint = rep(1, constants$numStates),
			
			p_star = rep(0.3, constants$numCounties),
			p_max = rep(0.35, constants$numStates),
			p_min_star = rep(0.25, constants$numStates)			
		)

		### monitors
		monitors <- names(inits)
		monitors <- monitors[!(monitors %in% c('z'))]
		monitors <- c(monitors, 'p_min')

	### model setup
	###############
	
		code <- nimble::nimbleCode({
			
			### likelihood
			for (i in 1:numCounties) {

				# detection
				y[i] ~ dbin(p[i], N[i])
				p[i] <- z[i] * p_star[i] + (1 - z[i]) * q
				p_star[i] ~ dunif(p_min[state[i]], p_max[state[i]])

				# occupancy
				z[i] ~ dbern(psi[i])
				psi[i] ~ dbeta(1, 1)
			
			}

			### priors
			q ~ dbeta(1, 10)
			
			# detection
			for (j in 1:numStates) {
				p_max[j] ~ dbeta(1, 1)
				p_min[j] <- p_min_star[j] * p_max[j]
				p_min_star[j] ~ dbeta(10, 1)
				qConstraint[j] ~ dconstraint(p_min[j] > q)
			}
			
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

		for (i in 1:constants$numStates) {

			node <- paste0('p_max[', i, ']')
			conf$removeSamplers(node)
			conf$addSampler(target=node, type='slice')
			
		}
				
		confBuild <- nimble::buildMCMC(conf)
		flush.console()
		compiled <- nimble::compileNimble(model, confBuild)
		flush.console()

	### MCMC
	########
		
		mcmc <- runMCMC(compiled$confBuild, niter = niter, nburnin = nburnin, nchains = nchains, inits = inits, progressBar = verbose, samplesAsCodaMCMC = TRUE, summary = FALSE, WAIC = TRUE)
		flush.console()

	### process model output
	########################

		mcmcModel <- .processBayesODM(shape=shape, mcmc=mcmc, effort=effort, detect=detect, hasNeighs=NULL)

		meta <- list(
			niter=niter,
			nburnin=nburnin,
			nchains=nchains,
			thin=thin,
			effort=effort,
			detect=detect,
			stateProv=stateProv,
			county=county,
			hasIslands=NA,
			na.rm=na.rm
		)
		
		mcmcModel <- c(mcmcModel, list(code=code), meta=list(meta))
		mcmcModel

}

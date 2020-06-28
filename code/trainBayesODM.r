#' Create a bias-corrected distribution model
#' 
#' This function implements a bayesLand model
#' @param x SpatialPolygonsDataFrame
#' @param effort Character, name of field in \code{x} indicating collection effort.
#' @param detect Character, name of field in \code{x} indicating number of detections of focal species.
#' @param stateProv Character, name of field with state/province.
#' @param county Character, name of field with county.
#' @param qGivenDetect Either \code{NULL} (default) or logical. If \code{TRUE} or \code{FALSE} then the model will try to estimate the rate of false detection (mistaken positive identifications). If \code{NULL}, then the model will not.  If \code{TRUE}, false detection rate \emph{q} reflects the probability of all collections in a location being of mistaken identification. If \code{FALSE} (default), then \emph{q} reflects the probability that for a county where the species does not truly occur, similar species (which could be mistaken for the focal species) do occur (and may or may not have been sampled and mis-identified).
#' @param covariate Either \code{NULL} (default) or a character or character vector. If not \code{NULL}, these are the name(s) field(s) with covariates.
#' @param niter Positive integer, number of MCMC iterations (including burn-in). The default is 2000, but this is often too low for most cases (i.e., try numbers in the range 10000, 100000, etc.).
#' @param nburnin Positive integer, number of burn-in samples (less than \code{niter}). The default is 1000, but this is often too low for most cases (i.e., try numbers in the range 10000, 100000, etc.).
#' @param nchains Positive integer, number of MCMC chains (default is 4).
#' @param thin Positive integer, number of MCMC samples by which to thin the results (default is 1; i.e., no thinning). To reduce memory requirements, you can use \code{thin} while increasing \code{niter} and/or \code{nburnin}.
#' @param na.rm Logical, if \code{TRUE} (default), then remove rows in \code{x} that have \code{NA} in any input field. If this is \code{FALSE} and \code{NA}s occur in any fields then an error may occur.
#' @param verbose Logical, if \code{TRUE} (default) display progress.
#' @param ... Arguments to pass to \code{\link[nimble]{configureMCMC}} and \code{\link[nimble]{runMCMC}}.
#' @return A list object with three elements, one with the MCMC chains, a second with the object \code{x} with model output appended to the data portion of the object, and a third with model information. The new columns in \code{x} represent:
#' \itemize{
#' 		\item{\code{cvScale}} (only if \code{covariate} is not \code{NULL}): Value of centered and scaled covariate
#' 		\item{\code{psi}} Probability of occurrence
#' 		\item{\code{psi95CI}} 95% credibility interval for probability of occurrence
#' 		\item{\code{p}} Probability of detecting the focal species assuming it is present
#' 		\item{\code{p95CI}} 95% credibility interval for probability of detection
#' }
#' @examples
#' @export

trainBayesODM <- function(
	x,
	effort,
	detect,
	stateProv,
	county,
	qGivenDetect = NULL,
	covariate = NULL,
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
	
		if (verbose) omnibus::say('Preparing data...')
		
		# remove missings
		if (na.rm) {
			ok <- complete.cases(x@data[ , c(effort, detect, stateProv, county, covariate)])
			x <- x[ok, ]
		}
			
		### CAR setup for counties
		##########################
		
			neighs <- spdep::poly2nb(x, queen=TRUE)
			
			hasNeighs <- !sapply(neighs, function(x) { length(x == 1) && (x == 0) })
			
			neighsList <- spdep::nb2WB(neighs)
			carAdjCounty <- neighsList$adj
			carWeightCounty <- neighsList$weights
			carNumCounty <- neighsList$num

		### inputs
		##########
		
			### data
			data <- list(
				y = x@data[ , detect],
				N = x@data[ , effort],
				carAdjCounty = carAdjCounty,
				carWeightCounty = carWeightCounty,
				carNumCounty = carNumCounty
			)
			
			# add covariates
			if (!is.null(covariate)) {
				for (i in seq_along(covariate)) {
					predName <- covariate[i]
					predVals <- as.numeric(scale(x@data[ , predName]))
					predList <- list(predVals)
					names(predList) <- paste0('predictor', i)
					data <- c(data, predList)
					x@data$DUMMY <- predVals
					names(x@data)[ncol(x@data)] <- paste0('predictor', i)
				}
			}

			### constants
			constants <- list(
				# falseDetect = as.integer((!is.null(qGivenDetect))),
				yConst = data$y,
				hasNeighs = as.numeric(hasNeighs),
				numCounties = nrow(x),
				numStates = length(unique(x@data[ , stateProv])),
				state = as.numeric(as.factor(x@data[ , stateProv])),
				lengthAdjCounties = length(carAdjCounty)
			)
			
			if (!is.null(qGivenDetect)) {
				if (qGivenDetect) {
					isOcc <- as.numeric(data$y > 0)
				} else {
					isOcc <- rep(1, nrow(x))
				}
				
				constants <- c(constants, list(isOcc=isOcc))
			}
		
			### initializations
			z <- as.numeric(data$y > 0)
			psi <- runif(constants$numCounties, 0.4, 0.6)
			p <- runif(constants$numCounties, 0.1, 0.3)
			
			inits <- list(
			
				psi = psi,
				p = p,
				logit_p = logit(p),
				
				p_stateMean = rnorm(constants$numStates, 0, 0.5),
				
				psi_island = rnorm(constants$numCounties, 0, 0.1),
				psi_islandMean = -0.1,
				
				psi_tau = 0.1,
				psi_car = rnorm(constants$numCounties)
				
			)

			if (!is.null(qGivenDetect)) {
				inits <- c(
					inits,
					list(
						z = z,
						q = 0.01
					)
				)
			}
			
			if (!is.null(covariate)) {
			
				betas <- list(
					psi_beta1 = 0.1,
					psi_beta1pow2 = 0.1,
					psi_beta2 = 0.1,
					psi_beta2pow2 = 0.1,
					psi_beta3 = 0.1,
					psi_beta3pow2 = 0.1
				)
				
				inits <- c(inits, betas)

			}
			
			### monitors
			monitors <- names(inits)
			monitors <- monitors[!(monitors %in% c('logit_p'))]
			if (is.null(qGivenDetect)) {
				monitors <- monitors[!(monitors %in% c('z'))]
			} else {
				monitors2 <- c('q', 'qAll')
			}
	
	### model setup
	###############
	
		if (verbose) omnibus::say('Model setup...')
		
		### model specification
		if (is.null(covariate)) {
			
			### WITHOUT covariate
			code <- nimble::nimbleCode({
				
				## likelihood
				for (i in 1:numCounties) {

					# detection
					y[i] ~ dbin(pstar[i], N[i])
					if (falseDetect == 1) {
						pstar[i] <- psi[i] * p[i] + (1 - psi[i]) * isOcc[i] * qAll[i] # viable but gives large q if isOcc is {0, 1}
						qAll[i] <- 1 - (1 - q)^yConst[i]
					} else {
						# pstar[i] <- z[i] * p[i] # viable but gives large q if isOcc is {0, 1}
						# z[i] ~ dbern(psi[i])
						pstar[i] <- psi[i] * p[i]
						z[i] ~ dbern(psi[i])
					}
					logit(p[i]) ~ dnorm(p_stateMean[state[i]], tau=0.001)
					

					# occupancy
					logit(psi[i]) <- hasNeighs[i] * psi_car[i] + (1 - hasNeighs[i]) * psi_island[i]
					psi_island[i] ~ dnorm(psi_islandMean, tau=0.001)
				
				}

				## priors

				# state effects
				for (j in 1:numStates) {
					p_stateMean[j] ~ dnorm(0, tau=0.001)
				}
				
				psi_islandMean ~ dnorm(0, tau=0.001)
				
				# false detection
				if (falseDetect == 1) q ~ dbeta(1, 2)
				
				# county occupancy CAR
				psi_tau ~ dgamma(0.001, 0.001)
				psi_car[1:numCounties] ~ dcar_normal(adj=carAdjCounty[1:lengthAdjCounties], weights=carWeightCounty[1:lengthAdjCounties], num=carNumCounty[1:numCounties], tau=psi_tau, c=3, zero_mean=0)
				
			})
			
		} else {
		
			### WITH covariate
			code <- nimble::nimbleCode({
				
				## likelihood
				for (i in 1:numCounties) {

					# detection
					y[i] ~ dbin(pstar[i], N[i])
					# if (falseDetect == 1) {
						# pstar[i] <- psi[i] * p[i] + (1 - psi[i]) * isOcc[i] * qAll[i] # viable but gives large q if isOcc is {0, 1}
						# qAll[i] <- 1 - (1 - q)^yConst[i]
					# } else {
						# pstar[i] <- z[i] * p[i] # viable but gives large q if isOcc is {0, 1}
						# z[i] ~ dbern(psi[i])
						pstar[i] <- psi[i] * p[i]
						# z[i] ~ dbern(psi[i])
					# }
					logit(p[i]) ~ dnorm(p_stateMean[state[i]], tau=0.001)
					

					# occupancy
					logit(psi[i]) <- hasNeighs[i] * psi_car[i] + (1 - hasNeighs[i]) * psi_island[i] + psi_beta1 * predictor1[i] + psi_beta1pow2 * (predictor1[i]^2) + psi_beta2 * predictor2[i] + psi_beta2pow2 * (predictor2[i]^2) + psi_beta3 * predictor3[i] * psi_beta3pow2 * (predictor3[i]^2)
					psi_island[i] ~ dnorm(psi_islandMean, tau=0.001)
				
				}

				## priors

				# state effects
				for (j in 1:numStates) {
					p_stateMean[j] ~ dnorm(0, tau=0.001)
				}
				
				psi_islandMean ~ dnorm(0, tau=0.001)
				
				# # false detection
				# if (falseDetect == 1) q ~ dbeta(1, 2)
				
				# psi covaraites
				psi_beta1 ~ dnorm(0, tau=0.001)
				psi_beta1pow2 ~ dnorm(0, tau=0.001)
				psi_beta2 ~ dnorm(0, tau=0.001)
				psi_beta2pow2 ~ dnorm(0, tau=0.001)
				psi_beta3 ~ dnorm(0, tau=0.001)
				psi_beta3pow2 ~ dnorm(0, tau=0.001)
				
				# county occupancy CAR
				psi_tau ~ dgamma(0.001, 0.001)
				psi_car[1:numCounties] ~ dcar_normal(adj=carAdjCounty[1:lengthAdjCounties], weights=carWeightCounty[1:lengthAdjCounties], num=carNumCounty[1:numCounties], tau=psi_tau, c=3, zero_mean=0)
				
			})
			
		} # if covariate
		
		
		### construct and compile model
		model <- nimble::nimbleModel(code=code, constants=constants, data=data, inits=inits, check=TRUE)
		flush.console()
		
		# conf <- nimble::configureMCMC(model, monitors = monitors, thin = thin, print = verbose, ...)
		conf <- nimble::configureMCMC(model, monitors = monitors, thin = thin, print = verbose)
		flush.console()
		
		## modify samplers
		conf$removeSamplers('psi_tau')
		conf$addSampler(target='psi_tau', type='slice')
		
		if (!is.null(qGivenDetect)) {
			conf$removeSamplers('q')
			conf$addSampler(target='q', type='slice')
		}

		if (!is.null(covariate)) {
			conf$removeSamplers('psi_beta')
			conf$addSampler(target='psi_beta', type='slice')
		}

		for (i in 1:constants$numCounties) {
			conf$removeSamplers(paste0('logit_p[', i, ']'))
			conf$addSampler(target=paste0('logit_p[', i, ']'), type='slice')
		
		}
		
		for (i in 1:constants$numCounties) {
			conf$removeSamplers(paste0('psi_island[', i, ']'))
			conf$addSampler(target=paste0('psi_island[', i, ']'), type='slice')
		
		}

		if (!is.null(covariate)) {
			for (i in 1:3) {
				node <- paste0('psi_beta', i)
				conf$removeSamplers(node)
				conf$addSampler(target=node, type='slice')

				node <- paste0('psi_beta', i, 'pow2')
				conf$removeSamplers(node)
				conf$addSampler(target=node, type='slice')
			}
		}
				
		confBuild <- nimble::buildMCMC(conf)
		flush.console()
		compiled <- nimble::compileNimble(model, confBuild)
		flush.console()

	### MCMC
	########
		
		if (verbose) omnibus::say('Modeling...')
		# mcmc <- runMCMC(compiled$confBuild, niter = niter, nburnin = nburnin, nchains = nchains, inits = inits, progressBar = verbose, samplesAsCodaMCMC = TRUE, summary = FALSE, ...)
		mcmc <- runMCMC(compiled$confBuild, niter = niter, nburnin = nburnin, nchains = nchains, inits = inits, progressBar = verbose, samplesAsCodaMCMC = TRUE, summary = FALSE)
		flush.console()

		# rm(model, compiled, confBuild)
		# gc()
		
	### process model output
	########################
	
		omnibus::say('Processing model output...')

		### remove unneeded stuff
		#########################

			mcmc <- list(samples = mcmc)

			# remove unneeded
			ignores <- c('logit_p[[]', 'z[[]')
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

						# conv <- rhatStats(mcmc$samples, rhatThresh=1.1, minConv=minUnconverged)
						# print(str(conv))
						# w <- as.Bwiqid(mcmc$samples)

						# caterplot(mcmc$samples, regex='p[[]')
						# caterplot(mcmc$samples, regex='psi[[]')

						# diagPlot(w, params='p')
						# diagPlot(w, params='psi')
						# diagPlot(w, params='q')
						# diagPlot(w, params='psi_tau')
			
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
			
			mcmc$summary$all.chains[ , 'Mean'] <- colMeans(samps)
			mcmc$summary$all.chains[ , 'Median'] <- apply(samps, 2, median)
			mcmc$summary$all.chains[ , 'St.Dev'] <- apply(samps, 2, sd)
			mcmc$summary$all.chains[ , '95%CI_low'] <- apply(samps, 2, quantile, probs=0.025)
			mcmc$summary$all.chains[ , '95%CI_upp'] <- apply(samps, 2, quantile, probs=0.975)
			mcmc$summary$all.chains[ , 'rhat'] <- wiqid::simpleRhat(mcmc$samples, n.chains=nchains)

			rm(samps); gc()

		### update shapefile
		####################

			psi <- mcmc$summary$all.chains[grepl(rownames(mcmc$summary$all.chains), pattern='psi[[]'), 'Mean']
			p <- mcmc$summary$all.chains[grepl(rownames(mcmc$summary$all.chains), pattern='p[[]'), 'Mean']
			
			uncer <- mcmc$summary$all.chains[ , '95%CI_upp'] - mcmc$summary$all.chains[ , '95%CI_low']
			names(uncer) <- rownames(mcmc$summary$all.chains)
			occUncer <- uncer[grepl(names(uncer), pattern='psi[[]')]
			pUncer <- uncer[grepl(names(uncer), pattern='p[[]')]

			x@data$psi <- psi
			x@data$occUncer <- occUncer
			x@data$p <- p
			x@data$pUncer <- pUncer
			
			names(x@data)[(ncol(x@data) - 3):ncol(x@data)] <- c('psi', 'psi95CI', 'p', 'p95CI')
			
		### collate
		###########
		
			meta <- list(
				niter=niter,
				nburnin=nburnin,
				nchains=nchains,
				thin=thin,
				effort=effort,
				detect=detect,
				stateProv=stateProv,
				county=county,
				covariate=covariate,
				hasIslands=any(!hasNeighs),
				qGivenDetect=qGivenDetect,
				na.rm=na.rm
			)
			
			attr(x, 'meta') <- meta
				
	### output
	list(mcmc=mcmc, code=code, shape=x, meta=meta)

}

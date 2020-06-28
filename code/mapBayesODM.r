#' Make display maps for a species
#'
#' Make display maps for a species.
#' @param species Name of species
#' @param family Family of species
#' @param focus2Sp Spatial data frame object
#' @param gadm1Sp Spatial object with states, does not have to be cropped to focus2Sp
#' @param mcmc MCMC list object from NIMBLE
#' @param minUnconverged Proportion of unconverged parameters that is acceptable for a model
#' @param outDir Directory to which to save image
#' @param isOk TRUE/FALSE if model is sufficiently converged. If \code{NULL} the this will be calculated.
#' @param ... Arguments to pass to text()
#' @return Write an image to disk.
#' @export

mapBayesODM <- function(
	species,
	family,
	focus2Sp,
	gadm1Sp,
	mcmc,
	minUnconverged,
	outDir,
	isOk = NULL,
	...
) {

	### formatting
	occCol <- 'forestgreen' # green--map color for occupancy
	occUncertCol <- 'red' # map color for uncertainty in occupancy
	knownOccCol <- 'steelblue3' # map color for counties with known occurrences
	effortCol <- '#6a51a3' # blue--map color number of collections
	
	lwd <- 0.4 # line width of political borders

	### assess convergence
	######################

		# if (TRUE) {
			# ignore <- c('z[[]', 'falseDetect[[]', 'constraintFalseDetect[[]')
			# conv <- rhatStats(mcmc$samples, rhatThresh = 1.1, minConv = minUnconverged, ignore = ignore)
		# } else {
			# if (FALSE) conv <- list(); conv$sufficient <- TRUE
		# }

	### prepare plot elements
	#########################

		occ2Sp <- focus2Sp[focus2Sp@data$detect > 0, ]
	
		# centroids of counties with occurrences
		focus2SpEa <- sp::spTransform(focus2Sp, enmSdm::getCRS('mollweide', TRUE))
		centsSpEa <- rgeos::gCentroid(focus2SpEa[focus2SpEa@data$detect > 0, ], byid=TRUE)
		centsSp <- sp::spTransform(centsSpEa, raster::projection(focus2Sp))
		
		# create plotting extent: focal region plus a buffer
		extEa <- raster::extent(focus2SpEa)
		extSpEa <- as(extEa, 'SpatialPolygons')
		raster::projection(extSpEa) <- enmSdm::getCRS('mollweide')
		
		extCentSpEa <- rgeos::gCentroid(extSpEa)
		maxDist_m <- rgeos::gDistance(extCentSpEa, extSpEa, hausdorff=TRUE)
		extSpEa <- rgeos::gBuffer(extSpEa, width=0.05 * maxDist_m)
		
		extSp <- sp::spTransform(extSpEa, raster::projection(focus2Sp))
		ext <- raster::extent(extSp)
		extSp <- as(ext, 'SpatialPolygons')
		projection(extSp) <- raster::projection(focus2Sp)

		# gadm1Sp <- gadm1Sp[gadm1Sp@data$NAME_1 %in% unique(focusSp@data$NAME_1), ]
		gadm1Sp <- sp::spTransform(gadm1Sp, sp::CRS(raster::projection(focusSp)))
		gadm1SpCrop <- raster::crop(gadm1Sp, extSp)

	### plot
	########

		if (is.null(isOk)) {
			conv <- rhatStats(mcmcModel$mcmc$samples, rhatThresh=1.1, minConv=minUnconverged)
			isOk <- conv$sufficient
		}

		ok <- if (isOk) { 'ok' } else { 'notOk'}
		filename <- paste0(tolower(family), '_', species)
		png(paste0(outDir, '/', filename, '.png'), width=1800, height=1100, res=300)
		
			# layout
			par(bg='white', oma=c(0.7, 0, 1.3, 0), mar=c(0, 2, 0, 0))
			lay <- matrix(c(1, 1, 2, 1, 1, 3), nrow=2, byrow=TRUE)
			layout(lay)

		### occupancy
		
			col <- occCol
			
			plot(gadm1SpCrop, col='gray90', border=NA)
			
			x <- focus2Sp@data$psi
			if (any(focus2Sp@data$detect == 0)) x[focus2Sp@data$detect > 0] <- NA
			xMaxNonDetect <- max(x, na.rm=TRUE)
			cols <- scales::alpha(col, x / xMaxNonDetect)
			plot(focus2Sp, col=cols, border='gray80', lwd=lwd, add=TRUE)
			plot(occ2Sp, col=knownOccCol, border='gray80', lwd=lwd, add=TRUE)
			plot(gadm1SpCrop, col=NA, border='gray60', lwd=lwd, add=TRUE)
			
			# legend
			labels <- seq(0, xMaxNonDetect, length.out=5)
			labels <- round(labels, 2)
			labels <- sprintf('%0.2f', labels)
			labels[1] <- '0'
			
			width <- 0.04
			height <- 0.9
			
			legTitle <- 'Occupancy\nprobability'
			legCex <- 0.62

			swatches <- list(
				list(swatchAdjY=c(0, 0.02), col=knownOccCol, border='gray', labels='Specimens')
			)
			
			legendary::legendGrad('left', inset=-0.01, title=legTitle, col=c('white', col), labels=labels, width=width, height=height, border='gray', boxBorder=NA, adjX=c(0.6, 1), adjY=c(0.06, 0.90), titleAdj=c(1.39, 0.95), labAdj=0.82, boxBg=NA, cex=legCex, lwd=0.5 * lwd, swatches=swatches, pos=2)
			
		### uncertainty in occupancy
		
			col <- occUncertCol
			
			plot(gadm1SpCrop, col='gray90', border=NA)
			
			x <- focus2Sp$psi95CI
			if (any(focus2Sp@data$detect == 0)) x[focus2Sp@data$detect > 0] <- NA
			xMaxNonDetect <- max(x, na.rm=TRUE)
			cols <- scales::alpha(col, x / xMaxNonDetect)
			plot(focus2Sp, col=cols, border='gray80', lwd=lwd, add=TRUE)
			plot(occ2Sp, col=knownOccCol, border='gray80', lwd=lwd, add=TRUE)
			plot(gadm1SpCrop, col=NA, border='gray60', lwd=lwd, add=TRUE)
			
			# legend
			labels <- seq(0, xMaxNonDetect, length.out=5)
			labels <- round(labels, 2)
			labels <- sprintf('%0.2f', labels)
			labels[1] <- '0'
			
			legCex <- 0.58
			width <- 0.12
			height <- 0.3
			
			legTitle <- 'Occupancy\nuncertainty\n(95% CI)'

			swatches <- list(
				list(swatchAdjY=c(0, 0.04), col=knownOccCol, border='gray', labels='Collections')
			)
			
			legendary::legendGrad('left', inset=-0.01, title=legTitle, col=c('white', col), labels=labels, width=width, height=0.925, border='gray', boxBorder=NA, adjX=c(0, 0.29), adjY=c(0.07, 0.81), titleAdj=c(0.55, 0.93), labAdj=0.2, boxBg=NA, cex=legCex, lwd=0.5 * lwd, swatches=swatches, pos=2)
			
		### effort
		
			title <- expression(paste('Total '*italic(family.)))
			legTitle <- 'Specimens'
			legCex <- 0.35
			col <- effortCol
			
			plot(gadm1SpCrop, col='gray90', border=NA)
			
			x <- focus2Sp@data$effort
			x <- log10(x + 1)
			cols <- scales::alpha(col, x / max(x))
			plot(focus2Sp, col=cols, border='gray80', lwd=lwd, add=TRUE)
			plot(gadm1SpCrop, col=NA, border='gray60', lwd=lwd, add=TRUE)
			points(centsSp, pch=16, cex=0.25)
			
			# legend
			labels <- seq(0, 1, by=0.2) * max(x)
			labels <- round(10^labels - 1, digits=1)
			labels <- sprintf('%.1f', labels)
			labels[1] <- '0'
			
			width <- 0.12
			height <- 0.3
			
			legCex <- 0.58
			width <- 0.12
			height <- 0.3
			
			legTitle <- paste('Total\n', genus)

			legendary::legendGrad('left', inset=-0.01, title=legTitle, col=c('white', col), labels=labels, width=width, height=0.925, border='gray', boxBorder=NA, adjX=c(0, 0.29), adjY=c(0, 0.84), titleAdj=c(0.55, 0.93), labAdj=0.2, boxBg=NA, cex=legCex, lwd=0.5 * lwd, pos=2)

		### titles
		
			main <- substitute(paste(family, ': ', italic(species)), env=list(species=species, family=family))
			mtext(text=main, at=c(0.01), outer=TRUE, cex=1, line=-0.3, adj=0)
			if (!isOk) mtext(text='Insufficient convergence', at=c(0.01), outer=TRUE, cex=1, line=-1.3, adj=0, col='red')
			mtext(text=date(), side=1, at=0.99, outer=TRUE, cex=0.20, line=-0.3, adj=1)
			
			if ('q' %in% rownames(mcmc$summary$all.chains)) {
			
				q <- round(mcmc$summary$all.chains['q', 'Mean'], 3)
				qLow <- round(mcmc$summary$all.chains['q', '95%CI_low'], 3)
				qHigh <- round(mcmc$summary$all.chains['q', '95%CI_upp'], 3)

				msg <- paste0('Probability of mistaken identification in any county with a single specimen: ', sprintf('%0.3f', q), ' (95% CI: ', sprintf('%0.3f', qLow), '-', sprintf('%0.3f', qHigh), ')')
				mtext(msg, side=1, at=0.01, outer=TRUE, cex=0.5, line=-0.3, adj=0)
			
			}
			
		dev.off()
	
}

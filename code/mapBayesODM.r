#' Make display maps for a species
#'
#' Make display maps for a species.
#' @param species Name of species
#' @param family Family of species
#' @param focus2Sp Spatial data frame object in equal-area projection
#' @param gadm1Sp Spatial object with states, does not have to be cropped to focus2Sp
#' @param q Vector of posterior samples of false positive rate of detection
#' @param waic WAIC object or \code{NULL}
#' @param loo LOO object or \code{NULL}
#' @param minUnconverged Proportion of unconverged parameters that is acceptable for a model
#' @param outDir Directory to which to save image
#' @param isOk TRUE/FALSE if model is sufficiently converged.
#' @param appendToName \code{NULL} or character to append to file name.
#' @param footer \code{NULL} or character to append to add to footer.
#' @param ... Arguments to pass to text()
#' @return Write an image to disk.
#' @export

mapBayesODM <- function(
	species,
	family,
	focus2Sp,
	gadm1Sp,
	q,
	waic,
	loo,
	minUnconverged,
	outDir,
	isOk,
	appendToName = NULL,
	footer = NULL,
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

		# centroids of counties with occurrences
		occ2Sp <- focus2Sp[focus2Sp@data$detect > 0, ]
		centsSp <- rgeos::gCentroid(focus2Sp[focus2Sp@data$detect > 0, ], byid=TRUE)
		
		# fix any potentially orphaned holes
		focus2Sp <- rgeos::gBuffer(focus2Sp, width=0, byid=TRUE)
		
		# create plotting extent: focal region plus a buffer
		ext <- raster::extent(focus2Sp)
		extSp <- as(ext, 'SpatialPolygons')
		raster::projection(extSp) <- raster::projection(focusSp)
		
		extCentSp <- rgeos::gCentroid(extSp)
		maxDist_m <- rgeos::gDistance(extCentSp, extSp, hausdorff=TRUE)
		extSp <- rgeos::gBuffer(extSp, width=0.05 * maxDist_m)
		
		ext <- raster::extent(extSp)
		extSp <- as(ext, 'SpatialPolygons')
		projection(extSp) <- raster::projection(focus2Sp)

		# crop GADM1
		extSp_inGadm <- sp::spTransform(extSp, sp::CRS(raster::projection(gadm1Sp)))
		gadm1SpCrop <- raster::crop(gadm1Sp, extSp_inGadm)
		gadm1SpCrop <- sp::spTransform(gadm1SpCrop, sp::CRS(raster::projection(focusSp)))

	### plot
	########

		# ok <- if (isOk) { 'ok' } else { 'notOk'}
		ok <- 'notAssessed'
		filename <- paste0(tolower(family), '_', species)
		if (!is.null(appendToName)) filename <- paste0(filename, '_', appendToName)
		png(paste0(outDir, '/', filename, '.png'), width=1800, height=1100, res=300)
		
			# layout
			par(bg='white', oma=c(0.7, 0.5, 1.3, 0), mar=c(0, 2, 0, 0))
			lay <- matrix(c(1, 1, 2, 1, 1, 3), nrow=2, byrow=TRUE)
			layout(lay)

		### occupancy
		
			col <- occCol
			
			plot(gadm1SpCrop, col='gray90', border=NA)
			
			x <- focus2Sp@data$psi
			if (any(focus2Sp@data$detect == 0)) x[focus2Sp@data$detect > 0] <- NA
			xMaxNonDetect <- max(x, na.rm=TRUE)
			cols <- scales::alpha(col, x / xMaxNonDetect)
			plot(focus2Sp, col='white', border=NA, add=TRUE)
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
			
			x <- focus2Sp$psi90CI
			if (any(focus2Sp@data$detect == 0)) x[focus2Sp@data$detect > 0] <- NA
			xMaxNonDetect <- max(x, na.rm=TRUE)
			cols <- scales::alpha(col, x / xMaxNonDetect)
			plot(focus2Sp, col='white', border=NA, add=TRUE)
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
			
			legTitle <- 'Occupancy\nuncertainty\n(90% CI)'

			swatches <- list(
				list(swatchAdjY=c(0, 0.04), col=knownOccCol, border='gray', labels='Specimens')
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
			plot(focus2Sp, col='white', border=NA, add=TRUE)
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
			
			legTitle <- paste('Total\n', family)

			legendary::legendGrad('left', inset=-0.01, title=legTitle, col=c('white', col), labels=labels, width=width, height=0.925, border='gray', boxBorder=NA, adjX=c(0, 0.29), adjY=c(0, 0.84), titleAdj=c(0.55, 0.93), labAdj=0.2, boxBg=NA, cex=legCex, lwd=0.5 * lwd, pos=2)

		### titles
		
			main <- substitute(paste(family, ': ', italic(species)), env=list(species=species, family=family))
			mtext(text=main, at=c(0.01), outer=TRUE, cex=1, line=-0.3, adj=0)
			if (!isOk) mtext(text='Insufficient convergence', at=c(0.01), outer=TRUE, cex=1, line=-1.3, adj=0, col='red')
			
			text <- paste(footer, ' | ', date())
			mtext(text=text, side=1, at=0.985, outer=TRUE, cex=0.35, line=-0.27, adj=1)
			
			qMean <- round(mean(q), 3)
			qLow <- round(quantile(q, 0.05), 3)
			qHigh <- round(quantile(q, 0.95), 3)

			msg1 <- paste0('WAIC: ', round(waic$estimates['waic', 1], 2), ' | LOO: ', round(loo$estimates['looic', 1], 2))
			msg2 <- paste0('Probability all specimens in a county in which ', species, ' is truely absent are mistaken identifications: ', sprintf('%0.3f', qMean), ' (90% CI: ', sprintf('%0.3f', qLow), '-', sprintf('%0.3f', qHigh), ')')
			mtext(msg1, side=1, at=0.005, outer=TRUE, cex=0.35, line=-0.8, adj=0)
			mtext(msg2, side=1, at=0.005, outer=TRUE, cex=0.35, line=-0.3, adj=0)
			
		dev.off()
	
}

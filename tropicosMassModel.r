### TROPICOS.Species Modeling
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org
### This script runs distribution models for each species in TROPICOS with sufficient data.
###
### source('C:/Ecology/Drive/R/tropicosMassModeling/tropicosMassModel.r')

### CONTENTS ###
### setup
### generalization
### constants
### get GADM ###
### retrieve species data from GBIF ###
### process species data from GBIF ###
### download species data from BIEN ###
### model ###

### model ###

#############
### setup ###
#############

	# set working directory
	setwd('C:/Ecology/Drive/R/tropicosMassModeling')
		
	# manage memory
	memory.limit(memory.limit() * 2^30)
	rm(list=ls())
	options(keep.source=FALSE) # manage memory
	gc()
	print('')
	print(date())

	options(stringsAsFactors=FALSE)
	
	these <- 1:20 # IP
	# these <- 101:200 # IP
	
	### settings
	############

		### model settings
		niter <- 20000 # burn-in + sampling
		nburnin <- 10000
		nchains <- 4
		thin <- 10

		### quality assurance settings
		minUnconverged <- 0.025 # minimum number of coefficients that can be uncoverged (at rhat > <-  1.1) for a model to be used
		minNumDupsToModel <- 3 # minimum number of duplicate records in a county necessary to generate a model
		minNumEffortsToModel <- 200 # minimum number of other species records that must be in the study region to model the focal species
		minPropCountiesWithEffortToModel <- 0.7 # minimum proportion of counties in the study region with any effort necessary to model the species
		minCountiesPerState <- 10 # minimum number of counties to include state in modeling if it has no detections
		expand <- 0.3 # to define study region create buffer that has a radius of largest distance from centroid to edge of occupied area times this amount
		minOccCounties <- 5 # minimum number of counties a species must appear in to model
		minNumRecords <- 20 # minimum number of records (across/within counties) to model

		saveNotOk <- TRUE # if TRUE then save results for models of species with insufficient convergence

	### use-specific generalization
	###############################
		
		longLat <- c('longitude', 'latitude')
		family <- 'Apocynaceae'
	
	
	# load/install packages
	if (!require(BIEN)) {
		install.packages('BIEN')
		library(BIEN)
	}
	
	if (!require(coda)) {
		install.packages('coda')
		library(coda)
	}
	
	if (!require(data.table)) {
		install.packages('data.table')
		library(data.table)
	}
	
	if (!require(mcmcplots)) {
		install.packages('mcmcplots')
		library(mcmcplots)
	}
	
	if (!require(nimble)) {
		install.packages('nimble')
		library(nimble)
	}
	
	if (!require(spdep)) {
		install.packages('spdep')
		library(spdep)
	}
	
	if (!require(wiqid)) {
		install.packages('wiqid')
		library(wiqid)
	}
	
	if (!require(stringi)) {
		install.packages('stringi')
		library(stringi)
	}
	
	if (!require(sp)) {
		install.packages('sp')
		library(sp)
	}
	
	if (!require(tidyverse)) {
		install.packages('tidyverse')
		library(tidyverse)
	}
	
	if (!require(terra)) {
		install.packages('terra')
		library(terra)
		# # if loading terra causes an error, use:
		# remotes::install_github('rspatial/terra')
		# library(terra)
	}
	
	
	if (!require(raster)) {
		install.packages('raster')
		library(raster)
	}
	
	if (!require(rgeos)) {
		install.packages('rgeos')
		library(rgeos)
	}
	
	if (!require(scales)) {
		install.packages('scales')
		library(scales)
	}
	
	if (!require(spocc)) {
		install.packages('spocc')
		library(spocc)
	}
	
	if (!require(taxize)) {
		install.packages('taxize')
		library(taxize)
	}
	
	if (!require(omnibus)) {
		remotes::install_github('adamlilith/omnibus')
		library(omnibus)
	}
	
	if (!require(legendary)) {
		remotes::install_github('adamlilith/legendary')
		library(legendary)
	}
	
	if (!require(enmSdm)) {
		remotes::install_github('adamlilith/enmSdm')
		library(enmSdm)
	}

	if (!require(birdsEye)) {
		remotes::install_github('adamlilith/birdsEye')
		library(birdsEye)
	}

	# get functions not part of package yet
	fxs <- listFiles('./code', pattern='.r')
	for (fx in fxs) source(fx)
	
	# load country ISO3 codes... this will eventually be in the "data" folder of the package
	load('./support/countryISO3.rda')
	
	
	omnibus::dirCreate('./data')
	omnibus::dirCreate('./regions')
	omnibus::dirCreate('./models')
		
omnibus::say('####################################################################')
omnibus::say('### get countries/states/counties shapefile and match to climate ###')
omnibus::say('####################################################################')

	# This section needs only done once. Results are saved to "regions" folder.

	# load world shapefiles
	if (file.exists('./regions/gadm2_worldclim.rda')) {
	
		load('./regions/gadm2_worldclim.rda')
		load('./regions/gadm1.rda')
		load('./regions/gadm0.rda')
		
	# retrieve GADM and match to environmental data
	} else {
		
		### get world level 2 (counties)
		
			iso3 <- raster::getData('ISO3')
			iso3 <- iso3$ISO3
			
			dir.create('./temp', showWarnings=FALSE)
			gadm2 <- raster::getData('GADM', country=iso3[1], level=2, path='./temp')
			for (country in iso3[2:length(iso3)]) {
				this <- tryCatch(
					expr={
						raster::getData('GADM', country=country, level=2, path='./temp')
					},
					error={
						function(x) FALSE
					}
				)
				if (class(this) != 'logical') gadm2 <- rbind(gadm2, this)
			}
			
		### get world level 1 (states/provinces)
			
			gadm1 <- raster::getData('GADM', country=iso3[1], level=1, path='./temp')
			for (country in iso3[2:length(iso3)]) {
				this <- tryCatch(
					expr={
						raster::getData('GADM', country=country, level=1, path='./temp')
					},
					error={
						function(x) FALSE
					}
				)
				if (class(this) != 'logical') gadm1 <- rbind(gadm1, this)
			}

		### get world level 0 (countries)
			
			gadm0 <- raster::getData('GADM', country=iso3[1], level=0, path='./temp')
			for (country in iso3[2:length(iso3)]) {
				this <- tryCatch(
					expr={
						raster::getData('GADM', country=country, level=0, path='./temp')
					},
					error={
						function(x) FALSE
					}
				)
				if (class(this) != 'logical') gadm0 <- rbind(gadm0, this)
			}

		### process GADM
		
			gadm2 <- gadm2[ , c('NAME_0', 'NAME_1', 'NAME_2')]
			gadm1 <- gadm1[ , c('NAME_0', 'NAME_1')]
			gadm0 <- gadm0[ , c('NAME_0')]
			
		### match GADM2 to environmental data
		
			# retrieve environmental data
			dir.create('./temp', showWarnings=FALSE)
			wc <- raster::getData('worldclim', var='bio', res=10, path='./temp')
			elevation_m <- raster::getData('worldclim', var='alt', res=10, path='./temp')
			names(elevation_m) <- 'elev_m'
			
			env <- raster::stack(elevation_m, wc)
			
			# fill NA cells near coasts to account for fact that some records may not fall in a cell near a coast
			name <- names(env)
			for (i in 1:nlayers(env)) {
				env[[i]] <- focal(env[[i]], w=matrix(1, nrow=3, ncol=3), fun=mean, na.rm=TRUE, NAonly=TRUE)
			}
			names(env) <- name
			
			# extract to shape
			envToGadm <- terra::extract(env, gadm2, weight=TRUE, normalizeWeights=TRUE)

			# match to shape taking area-weighted mean value
			gadm2@data$elevation_m <- NA
			for (i in seq_along(envToGadm)) {
				x <- envToGadm[[i]][ , 'elev_m']
				wgt <- envToGadm[[i]][ , 'weight']
				gadm2@data$elevation_m[i] <- sum(x * wgt, na.rm=TRUE)
			}

			for (bio in 1:19) {

				thisBio <- paste0('bio', bio)
				gadm2@data$DUMMY <- NA
				colnames(gadm2@data)[ncol(gadm2@data)] <- thisBio
				for (i in seq_along(envToGadm)) {
					x <- envToGadm[[i]][ , thisBio]
					wgt <- envToGadm[[i]][ , 'weight']
					gadm2@data[i, thisBio] <- sum(x * wgt, na.rm=TRUE)
				}
				
			}
	
		### calculate area of states and counties
		
			gadm2Ea <- sp::spTransform(gadm2, enmSdm::getCRS('mollweide', TRUE))
			gadm2@data$areaKm2 <- rgeos::gArea(gadm2Ea, byid=TRUE) / 1000^2
		
			gadm1Ea <- sp::spTransform(gadm1, enmSdm::getCRS('mollweide', TRUE))
			gadm1@data$areaKm2 <- rgeos::gArea(gadm1Ea, byid=TRUE) / 1000^2
		
		### clean up
			
			files <- listFiles('./temp', recursive=TRUE)
			for (file in files) file.remove(file)
			unlink('./temp', recursive=TRUE)

		### save
		
			save(gadm2, file='./regions/gadm2_worldclim.rda')
			save(gadm1, file='./regions/gadm1.rda')
			save(gadm0, file='./regions/gadm0.rda')
	
			rm(wc, env, elevation_m, area_km2, envToGadm, gadm2Ea)
			gc()
	
	}
	
# omnibus::say('#######################################')
# omnibus::say('### retrieve species data from GBIF ###')
# omnibus::say('#######################################')
	
	# # This section would ideally load data from TROPICOS using an API and from GBIF.
	# # In lieu of that, we are just going to download GBIF data.
	# # At present, this code does not usually work because of a bug where R confuses URLs in cells as column names.
	# # The only way I know of to manage this is to do a manual download.
	
	# ### get genera and species in this family
	# omnibus::say('Retreiving species names in this family...')
	
	# famInfo <- taxize::get_gbifid_(sciname=family, method='backbone')
	# famInfo <- famInfo[[1]]
	# if (famInfo$status != 'ACCEPTED' | famInfo$matchtype != 'EXACT' | famInfo$kingdom != 'Plantae') stop('This family is not accepted by GBIF.')
	
	# genera <- taxize::children(family, db='itis')
	# genera <- as.data.frame(genera[[1]])
	
	# species <- taxize::children(genera$tsn[1], db='itis')
	# speciesList <- species$taxonname
	# if (nrow(genera) > 1) {
		# for (i in 2:nrow(genera)) {
			# thisSpecies <- taxize::children(genera$tsn[i], db='itis')
			# thisSpecies <- as.data.frame(thisSpecies[[1]])
			# speciesList <- c(speciesList, thisSpecies$taxonname)
		# }
	# }
	
	# gbifFields <- c('key', 'family', 'species', 'year', 'month', 'day', 'country', 'stateProvince', 'county', 'verbatimLocality', 'locality', 'habitat', 'decimalLongitude', 'decimalLatitude', 'coordinateUncertaintyInMeters')
		
	# now <- substr(as.character(Sys.time()), 1, 4)
	
	# occs <- data.frame()
	# for (thisSpecies in speciesList) {
	
		# omnibus::say(thisSpecies)
	
		# thisOccs <- rgbif::occ_data(
			# scientificName = thisSpecies,
			# year=paste0('1950,', now),
			# limit=10,
			# fields=gbifFields
		# )
			
		# # thisOccs <- spocc::occ(
			# # query=thisSpecies,
			# # from='gbif',
			# # # limit=1000000,
			# # limit=10,
			# # date=c('1950-01-01', now),
			# # gbifopts=list(
				# # fields=c('species', 'year')
			# # )
				
		# # )
		
		# thisOccs <- as.data.frame(thisOccs$data)
		# thisOccs <- thisOccs[ , gbifFields]
		# occs <- rbind(occs, thisOccs)
		
	# }

	# ### remove suspicious records
	# occs$use <- TRUE
	
	# verbotens <- c('garden', 'experiment', 'greenhouse', 'hothouse', 'jardin', 'purchase', 'bought', 'wedding', 'arboret')
	# hab <- tolower(occs$habitat)
	# loc <- tolower(occs$locality)
	# verbLoc <- tolower(occs$verbatimLocality)
	# bads <- integer()
	# for (verboten in verbotens) {
		# bad <- which(grepl(hab, pattern=verboten) | grepl(loc, pattern=verboten) | grepl(verbLoc, pattern=verboten))
		# bads <- c(bads, bad)
	# }
	
	# if (length(bads) > 0) {
		# bads <- sort(unique(bads))
		# occs$use[bads] <- FALSE
	# }

	# omnibus::dirCreate('./data')
	# save(occs, file=paste0('./data/', tolower(family), '_occurrences.rda'))

# omnibus::say('######################################')
# omnibus::say('### process species data from GBIF ###')
# omnibus::say('######################################')

	# omnibus::say('Loading manual download from GBIF. Doing this instead of automatic download for now because of issue with URLs in column names.', breaks=80)
	
	# if (file.exists(paste0('./data/', tolower(family), '_occurrencesGadm.rda')) & file.exists(paste0('./data/', tolower(family), '_occurrences.rda'))) {
	
		# load(paste0('./data/', tolower(family), '_occurrencesGadm.rda'))
		# load(paste0('./data/', tolower(family), '_occurrences.rda'))
	
	# } else {
		
		# omnibus::say('Locating records to Earth...')
		
		# occs <- data.table::fread('./data/apocynaceae_northAmerica/occurrence.txt', header=TRUE)
		# occs <- as.data.frame(occs)
		
		# occs$newState <- NA
		# occs$newCounty <- NA

		# ## occurrences with coordinates
		
			# occsWithCoords <- which((!is.na(occs$decimalLongitude) | !is.na(occs$decimalLatitude)) & occs$coordinateUncertaintyInMeters %<na% 10000)

			# if (length(occsWithCoords) > 0) {
				
				# ext <- terra::extract(gadm2, occs[occsWithCoords, longLat])
				# occs$newState[occsWithCoords] <- ext$NAME_1
				# occs$newCounty[occsWithCoords] <- ext$NAME_2
				
			# }
				
		# # occurrences with just state/county names
		
			# occsNoCoords <- which(is.na(occs$decimalLongitude) | !is.na(occs$decimalLatitude))
		
			# if (length(occsNoCoords) > 0) {
				
				# occs$newState[occsNoCoords] <- occs$stateProvince[occsNoCoords]
				# occs$newCounty[occsNoCoords] <- occs$county[occsNoCoords]
				
				# occs$newState[occsNoCoords] <- stringi::stri_trans_general(occs$newState[occsNoCoords], 'latin-ascii')
				# occs$newCounty[occsNoCoords] <- stringi::stri_trans_general(occs$newCounty[occsNoCoords], 'latin-ascii')
				
			# }
			
		# ### remove unusable records

		# if (any(occs$species == '')) occs <- occs[occs$species != '', ]
		# if (any(occs$newState == '')) occs <- occs[occs$newState != '', ]
		# if (any(occs$newCounty == '')) occs <- occs[occs$newCounty != '', ]
		# if (any(is.na(occs$species))) occs <- occs[!is.na(occs$species), ]
		# if (any(is.na(occs$newState))) occs <- occs[!is.na(occs$newState), ]
		# if (any(is.na(occs$newCounty))) occs <- occs[!is.na(occs$newCounty), ]
		
		# if (any(occs$year < 1950)) occs <- occs[occs$year >= 1950, ]
		
		# ### remove cultivated/unnatural specimens
		# verbotens <- c('garden', 'experiment', 'greenhouse', 'hothouse', 'jardin', 'purchase', 'bought', 'wedding', 'arboret')
		# hab <- tolower(occs$habitat)
		# loc <- tolower(occs$locality)
		# verbLoc <- tolower(occs$verbatimLocality)
		# bads <- integer()
		# for (verboten in verbotens) {
			# bad <- which(grepl(hab, pattern=verboten) | grepl(loc, pattern=verboten) | grepl(verbLoc, pattern=verboten))
			# bads <- c(bads, bad)
		# }
		
		# if (length(bads) > 0) {
			# bads <- sort(unique(bads))
			# occs <- occs[-bads, ]
		# }

		# if (nrow(occs) > 0) {
			
			# # occs <- occs[ , gbifFields, drop=FALSE]
			# occs <- occs[ , c('gbifID', 'family', 'species', 'year', 'month', 'day', 'countryCode', 'stateProvince', 'county', 'verbatimLocality', 'locality', 'habitat', 'decimalLongitude', 'decimalLatitude', 'coordinateUncertaintyInMeters', 'newState', 'newCounty')]

			# occsGadm2 <- matchEffortToGadm(shape=gadm2, df=occs, dfStateProvField='newState', dfCountyField='newCounty')
			
			# save(occs, file=paste0('./data/', tolower(family), '_occurrences.rda'))
			# save(occsGadm2, file=paste0('./data/', tolower(family), '_occurrencesGadm.rda'))
			
		# }
			
	# }

omnibus::say('#######################################')
omnibus::say('### download species data from BIEN ###')
omnibus::say('#######################################')

	if (file.exists(paste0('./data/', tolower(family), '_occurrencesGadm.rda'))  & file.exists(paste0('./data/', tolower(family), '_occurrences.rda'))) {
	
		load(paste0('./data/', tolower(family), '_occurrences.rda'))
		load(paste0('./data/', tolower(family), '_occurrencesGadm.rda'))

	} else {

		occs <- BIEN::BIEN_occurrence_family(
		  family = family,
		  cultivated = FALSE,
		  only.new.world = TRUE,
		  observation.type = FALSE,
		  all.taxonomy = FALSE,
		  native.status = FALSE,
		  natives.only = TRUE,
		  political.boundaries = FALSE,
		  collection.info = FALSE
		)

		maxRows <- 10000
		sets <- ceiling(nrow(occs) / maxRows)
		
		matches <- data.frame(
			ntate = rep(NA, nrow(occs)),
			county = rep(NA, nrow(occs))
		)
		
		minimalGadm2 <- gadm2
		minimalGadm2 <- minimalGadm2[ , c('NAME_1', 'NAME_2')]
		
		for (set in 1:sets) {
		
			say('Extracting set ', set, ' of ', sets, ' of records.')
		
			rows <- ((set - 1) * maxRows + 1):((set - 1) * maxRows + maxRows)
			theseOccs <- occs[rows, longLat]
			ext <- terra::extract(minimalGadm2, theseOccs)
		
			matches$state[rows] <- ext$NAME_1
			matches$county[rows] <- ext$NAME_2
		
		}
		
		rm(minimalGadm2); gc()

		occs <- rbind(occs, matches)
		
		occsGadm2 <- matchEffortToGadm(shape=gadm2, df=occs, dfStateProvField='state', dfCountyField='county')
		
		save(occs, file=paste0('./data/', tolower(family), '_occurrences.rda'))
		save(occsGadm2, file=paste0('./data/', tolower(family), '_occurrencesGadm.rda'))

	}
	
omnibus::say('#############')
omnibus::say('### model ###')
omnibus::say('#############')

	# list of species
	speciesList <- sort(unique(occs$scrubbed_species_binomial))
	if (any(speciesList == '')) speciesList <- speciesList[-which(speciesList == '')]
	if (any(is.na(speciesList))) speciesList <- speciesList[-which(is.na(speciesList))]

	# gadm2 <- gadm2[gadm2@data$NAME_0 %in% focalCountries, ] # OPTIONAL!
	# gadm1 <- gadm2[gadm1@data$NAME_0 %in% focalCountries, ] # OPTIONAL!
	# gadm0 <- gadm2[gadm0@data$NAME_0 %in% focalCountries, ] # OPTIONAL!
	
	omnibus::dirCreate(paste0('./models/ok'))
	omnibus::dirCreate(paste0('./models/notOk'))

	### model each species
	# for (species in speciesList) {
	for (species in speciesList[these]) {
	# for (species in 'Ascelpias albicans') {
	# for (species in 'Vallesia glabra') {
	# for (species in 'Amsonia illustris') {

		omnibus::say(species, ' ', date(), level=2)
	
		test1 <- sum(occs$scrubbed_species_binomial == species, na.rm=TRUE) >= minNumRecords
	
		if (test1) {
		
			say('Insufficient records for modeling.')
			
		} else {

			# collate species into GADM
			occsGadm2 <- matchDetectToGadm(shape=occsGadm2, df=occs, dfStateProvField='state', dfCountyField='county', dfSpeciesField='scrubbed_species_binomial', species=species)

			numOccCounties <- sum(occsGadm2@data$detect > 0)
			
			test1 <- numOccCounties >= minOccCounties
			test2 <- any(occsGadm2@data$detect >= minNumDupsToModel)
			# test3 <- any(occsGadm2@data$NAME_0[occsGadm2@data$detect > 0] %in% focalCountries) # OPTIONAL!

			### model species
			if (!test1 | !test2) {
			
				say('Insufficient detections for modeling.')
			
			} else {

				famSpp <- paste0(tolower(family), '_', tolower(gsub(species, pattern=' ', replacement='_')))
		
				### study region definition
				###########################
				
					omnibus::say('   Defining modeling region...')

					focusSp <- occsGadm2[ , c('NAME_0', 'NAME_1', 'NAME_2', 'effort', 'detect', paste0('bio', 1:19), 'elevation_m')]
				
					focusSp <- getGeogFocus(
						x = focusSp,
						detect = 'detect',
						expand = expand,
						upper = 'NAME_1',
						minLowerPolysPerUpper = minCountiesPerState
					)
					
					# if there is sufficient effort
					test1 <- sum(focusSp@data$effort, na.rm=TRUE) > minNumEffortsToModel
					test2 <- sum(focusSp@data$effort > 0, na.rm=TRUE) / nrow(focusSp) >= minPropCountiesWithEffortToModel
					
					if (!test1 | !test2) {
					
						omnibus::say('Insufficient effort to model this species.')
						
					} else {
						
						# make an equal-area projection for this region
						
							ext <- raster::extent(focusSp)
							lat0 <- mean(c(ext@ymin, ext@ymax))
							lat1 <- ext@ymin
							lat2 <- ext@ymax
							long <- mean(c(ext@xmin, ext@xmax))
							
							eaProj <- paste0('+proj=lcc +lat_1=', lat1, ' +lat_2=', lat2, ' +lat_0=', lat0, ' +lon_0=', long, ' +x_0=0 +y_0=0 +ellps=GRS80 +datum=WGS84 +units=m +no_defs')
							
							focusSp <- sp::spTransform(focusSp, sp::CRS(eaProj)) # OPTIONAL!
						
					### PCA
					#######
					
						pca <- prcomp(focusSp@data[ , c(paste0('bio', 1:19), 'elevation_m')], center=TRUE, scale=TRUE)
						focusSp@data$pc1 <- pca$x[ , 1]
						focusSp@data$pc2 <- pca$x[ , 2]
						focusSp@data$pc3 <- pca$x[ , 3]
						
					### model
					#########
					
						mcmcModel <- trainBayesODM_simple(
							x=focusSp,
							effort='effort',
							detect='detect',
							stateProv='NAME_1',
							county='NAME_2',
							covariate=c('pc1', 'pc2', 'pc3'),
							qGivenDetect=NULL,
							niter=niter,
							nburnin=nburnin,
							nchains=nchains,
							thin=thin,
							na.rm=TRUE,
							verbose=TRUE
						)
						
					### process model output
					########################
					
						conv <- rhatStats(mcmcModel$mcmc$samples, rhatThresh=1.1, minConv=minUnconverged)
						isOk <- conv$sufficient
						ok <- if (isOk) { 'ok' } else { 'notOk' }
						
						print(str(conv))
						flush.console()

						### output
						##########
						
							if (saveNotOk | isOk) {
							
							outDir <- paste0('./models/', ok, '/', famSpp)
							omnibus::dirCreate(outDir)
							
							# save(mcmcModel, file=paste0(outDir, '/', famSpp, '_mcmcModel.rda'))

							### stats
							#########
				
								say('Statistics...')
							
								detecteds <- focusSp@data$detect
								efforts <- focusSp@data$effort
								
								sink(paste0(paste0(outDir, '/', famSpp, '.txt')))
								
									cat(species, '\n')
									cat(date(), '\n')
									cat('enmSdmBayes model. For more information,\nvisit www.earthSkySea.org/enmSdmBayes.\n\n')
									
									cat('INPUT\n')
									cat('number of occupied counties ', numOccCounties, '\n')
									cat('number of counties in focal region ', nrow(focusSp), '\n')
									cat('number of detections ', sum(detecteds), '\n')
									cat('maximum number of detections in a county ', max(detecteds), '\n')
									cat('number of effort collections ', sum(efforts), '\n')
									cat('average detected:effort ratio ', round(mean(detecteds / efforts, na.rm=TRUE), 4), '\n\n')
									
									cat('MODEL SETTINGS\n')
									cat('niter ', niter, '\n')
									cat('nburnin ', nburnin, '\n')
									cat('nchains ', nchains, '\n')
									cat('thin ', thin, '\n')
									cat('total samples ', (niter - nburnin) / thin, ' per chain\n\n')

									cat('MODEL OUTPUT\n')
									if (exists('WAIC', where=mcmcModel$mcmc)) cat('WAIC ', round(mcmcModel$mcmc$WAIC, 3), '\n')
									
									cat('Proportion of non-converged parameters (rhat > 1.1) ', conv$propUnconv, '\n\n')

									cat('Summary of rhats:', '\n\n')
									print(summary(conv$rhat))
									
									cat('\nUnconverged coefficients (rhat > 1.1):', '\n')
									print(conv$unconv)

									cat('\nNA coefficients (rhat > 1.1):', '\n')
									print(conv$nas)

								sink()
								
							### diagnostics plots
							#####################
							
								say('Diagnostic plots...')
							
								occs2Sp <- focusSp[focusSp@data$detect > 0, ]
								centsSp <- gCentroid(occs2Sp, byid=TRUE)
						
								png(paste0(outDir, '/', famSpp, '_diagnostics.png'), width=1600, height=1200, res=120)
								
									par(mfrow=c(2, 5), bg='gray90', oma=c(2, 0, 2, 0), mar=c(3, 4, 4, 2), cex.main=1.2)

									# effort
									scale <- rep(NA, nrow(focusSp))
									scale <- log10(focusSp@data$effort + 1) / max(log10(focusSp$effort + 1), na.rm=TRUE)
									cols <- scales::alpha('red', scale)
									cols[is.na(focusSp@data$effort)] <- 'gray'
									plot(focusSp, main='effort (log)', col=cols, border=NA)
									points(centsSp, pch=16, cex=0.2)

									# detections
									scale <- rep(NA, nrow(focusSp))
									scale <- log10(focusSp@data$detect + 1) / max(log10(focusSp@data$detect + 1), na.rm=TRUE)
									cols <- scales::alpha('black', scale)
									cols[is.na(focusSp$detect)] <- 'gray'
									plot(focusSp, main='detections', col=cols, border='gray80')

									# detectability
									thisIndex <- which(grepl(rownames(mcmcModel$mcmc$summary$all.chains), pattern='p[[]'))
									out <- mcmcModel$mcmc$summary$all.chains[thisIndex, 'Mean']
									cols <- ifelse(is.na(out), 'gray', scales::alpha('blue', out))
									plot(focusSp, main='detectability', col=cols, border=NA)
									points(centsSp, pch=16, cex=0.2)
									
									# psi
									thisIndex <- which(grepl(rownames(mcmcModel$mcmc$summary$all.chains), pattern='psi[[]'))
									out <- mcmcModel$mcmc$summary$all.chains[thisIndex, 'Mean']
									cols <- ifelse(is.na(out), 'gray', scales::alpha('darkgreen', out))
									plot(focusSp, main='occupancy', col=cols, border=NA)
									points(centsSp, pch=16, cex=0.2)

									states <- unique(focusSp@data$NAME_1)
									caterplot(mcmcModel$mcmc$samples, 'p_stateMean', reorder=FALSE, labels=states); title(main='state p')
									mcmcplots::caterplot(mcmcModel$mcmc$samples, regex='p[[]', reorder=TRUE); title(main='p')
									mcmcplots::caterplot(mcmcModel$mcmc$samples, regex='psi[[]', reorder=TRUE); title(main='psi')
									mcmcplots::caterplot(mcmcModel$mcmc$samples, regex='psi_tau'); title(main='psi tau')
									mcmcplots::caterplot(mcmcModel$mcmc$samples, regex='psi_car[[]', reorder=FALSE); title(main='psi CAR')
									mcmcplots::caterplot(mcmcModel$mcmc$samples, regex='psi_beta', reorder=FALSE); title(main='psi betas')

									if (isOk) {
										main <- paste0(family, ': ', species)
										col <- 'black'
									} else {
										main <- paste0(family, ': ', species, ' UNCONVERGED (rhat >1.1): ', sprintf('%0.1f', round(100 * conv$propUnconv, 1)), '%')
										col <- 'red'
									}
									
									title(main=main, cex.main=1.8, outer=TRUE, col=col)
									title(sub=date(), cex.sub=0.7, outer=TRUE, line=-0.3)
									
								dev.off()

							### save shapefile
							##################

								say('Shapefile...')
							
								out2Sp <- mcmcModel$shape
								out2Sp@data <- out2Sp@data[ , c('NAME_0', 'NAME_1', 'NAME_2', 'effort', 'detect', 'psi', 'psi95CI', 'p', 'p95CI')]
								names(out2Sp@data) <- c('country', 'state', 'county', 'effort', 'detect', 'psi', 'psi95CI', 'p', 'p95CI')

								if (isOk) {
									
									dirCreate(outDir, '/shapefiles')
									
									raster::shapefile(out2Sp, paste0(outDir, '/shapefiles/', famSpp), overwrite=TRUE)

									sink(paste0(outDir, '/shapefiles/', famSpp, '_README.txt'))
									
										say(family, ': ', species)
										say(date(), post=2)
										
										say('This shapefile contains the output of enmSdmBayesLand, a bias-corrected\nBayesian distribution model for ', species, '. The model is\nbased on an occupancy-detection framework, so attempts to correct for bias\nin collection (e.g., too much/little collection effort in certain locations).\nThe model uses non-duplicate collections of ', family, ' as a proxy of\ncollection effort and the number of non-duplicate collections of \n', species, ' as an index of "detectability."', post=2)
										
										say('The "sampling unit" in the model is a county or equivalent--the model\noutput is the probability of detecting the species in a county given that\nit is present in the county, and the probability of occurrence in the\ncounty. enmSdmBayesLand accounts for differences in detectability between\nstates/provinces, area of counties, and spatial autocorrelation in\noccurrence. It also provides estimates of uncertainty in detectability\nand occurrence.', post=2)
										
										say('The shapefile is based on GADM Version 3.6 (https://gadm.org/) and contains\nthe following fields:')
										say('* country: Country')
										say('* state: State/province')
										say('* county: County, parish, etc.')
										say('* areaKm2: Area of the county in km2')
										say('* effort: Number of non-duplicate records of ', family)
										say('* detect: Number of non-duplicate collections of ', species)
										say('* psi: Probability of occurrence of the species')
										say('* psi95CI: Width of the 95% credibility intrerval for psi')
										say('* p: Probability of collection of ', species, ' *given that it\nis present*')
										say('* p95CI: Width of the 95% credibility interval for p')
										
										say('enmSdmBayesLand was created by Camilo Sanín and Adam B. Smith.', pre=1)
										say('This shapefile and any model output are provided under the\nGNU GPLv3 license.')
										
									sink()						

									toZip <- listFiles(paste0(outDir, '/shapefiles'), pattern=famSpp)
									a <- zip::zipr(paste0(outDir, '/', famSpp, '.zip'), toZip)
									
									done <- file.remove(paste0(outDir, '/shapefiles/', famSpp, '.dbf'))
									done <- file.remove(paste0(outDir, '/shapefiles/', famSpp, '.prj'))
									done <- file.remove(paste0(outDir, '/shapefiles/', famSpp, '.shp'))
									done <- file.remove(paste0(outDir, '/shapefiles/', famSpp, '.shx'))
									done <- file.remove(paste0(outDir, '/shapefiles/', famSpp, '_README.txt'))
									unlink(paste0(outDir, '/shapefiles'), recursive=TRUE)
									
								} # if sufficient convergence
									
						### display plots
						#################

							say('Maps...')
						
							mapBayesODM(
								species=species,
								family=family,
								focus2Sp=out2Sp,
								gadm1Sp=gadm1,
								mcmc=mcmcModel$mcmc,
								minUnconverged=minUnconverged,
								isOk=isOk,
								outDir=outDir#,
								# pos=4
							)
							
					} # if model is sufficient or saving output from insufficient models

					### clean up
					
						rm(mcmcModel, focus2Sp, out2Sp)
						gc()
						
				} # if there is sufficient effort
						
			} # if species has sufficient detections
			
		} # if species has sufficient number of records
			
	} # next species

say('DONE!', deco='%', level=1)

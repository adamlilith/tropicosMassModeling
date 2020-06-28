#' Collate number of detections and effort into a GADM SpatialPolygonsDataFrame
#'
#' This function merges information on detections and effort from a data frame into a SpatialPolygonsDataFrame object from GADM.
#' @param shape SpatialPolygonsDataFrame in GADM level 2 format (www.gadm.org).
#' @param df Data frame
#' @param stateProv Name of field in \code{df} with stateProv/province names.
#' @param county Name of field in \code{df} with county names.
#' @param species Name of field in \code{df} with species names.
#' @param speciesName Name of focal species or \code{NULL} (do not match a particular species... just total collections in a polygon).
#' @examples

matchDfToGadm <- function(
	shape,
	df,
	stateProv = 'NAME_1',
	county = 'NAME_2',
	species = 'species',
	speciesName = NULL
) {

	if (!is.null(speciesName)) {
		if (!(speciesName %in% df[ , species])) warning('The species does not appear in this data frame.')
		shape@data$detect <- 0
		unfoundSpecies <- 0 # number of focal species' records that cannot be matched
	}

	shape@data$effort <- 0
	unfoundAll <- 0 # number of all records that cannot be matched

	# lower case everything
	dfStateProv <- tolower(df[ , stateProv])
	dfCounty <- tolower(df[ , county])
	
	shapeStateProv <- tolower(shape@data$NAME_1)
	shapeCounty <- tolower(shape@data$NAME_2)
	
	# replace characters with diacritics... makes matching easier
	shapeStateProv <- stringi::stri_trans_general(shapeStateProv, 'latin-ascii')
	shapeCounty <- stringi::stri_trans_general(shapeCounty, 'latin-ascii')
	
	dfStateProv <- stringi::stri_trans_general(dfStateProv, 'latin-ascii')
	dfCounty <- stringi::stri_trans_general(dfCounty, 'latin-ascii')

	# bad strings
	dfCounty <- gsub(dfCounty, pattern=' county', replacement='')
	dfCounty <- gsub(dfCounty, pattern=' cty.', replacement='')
	dfCounty <- gsub(dfCounty, pattern=' cty', replacement='')
	dfCounty <- gsub(dfCounty, pattern=' city', replacement='')
	dfCounty <- gsub(dfCounty, pattern='city of ', replacement='')
	dfCounty <- gsub(dfCounty, pattern='st. ', replacement='saint ')
	dfCounty <- gsub(dfCounty, pattern='ste. ', replacement='saint ')
	
	dfStateProv <- gsub(dfStateProv, pattern=' ', replacement='')
	dfCounty <- gsub(dfCounty, pattern=' ', replacement='')
	
	shapeStateProv <- gsub(shapeStateProv, pattern=' ', replacement='')
	shapeCounty <- gsub(shapeCounty, pattern=' ', replacement='')
	
	for (i in 1:nrow(df)) {
	
		# match data frame state/province and county to shape's state/province and county
		shapeIndex <- which(shapeCounty == dfCounty[i] & shapeStateProv == dfStateProv[i])
		
		# no geographic match
		if (length(shapeIndex) == 0) {
		
			unfoundAll <- unfoundAll + 1
			if (!is.null(speciesName)) {
				if (!is.na(df[i, species]) && df[i, species] == speciesName) {
					unfoundSpecies <- unfoundSpecies + 1
				}
			}
	
		# geographic match and is the focal species
		} else  {
		
			shape@data$effort[shapeIndex] <- shape@data$effort[shapeIndex] + 1

			if (!is.null(speciesName)) {
				if (!is.na(df[i, species]) && df[i, species] == speciesName) {
					shape@data$detect[shapeIndex] <- shape@data$detect[shapeIndex] + 1
				}
			}
			
		}
		
	}
	
	if (unfoundAll > 0) warning(unfoundAll, ' records were not matched in the shape.')
	if (!is.null(speciesName) && unfoundSpecies > 0) warning(unfoundSpecies, ' records of ', speciesName, ' were not matched in the shape.')

	shape

}

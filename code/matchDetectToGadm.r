#' Collate number of detections into a GADM SpatialPolygonsDataFrame
#'
#' This function merges information on detections of a focal species from a data frame into a SpatialPolygonsDataFrame object from GADM.
#' @param shape SpatialPolygonsDataFrame in GADM level 2 format (www.gadm.org).
#' @param df Data frame
#' @param species Name of species.
#' @param dfStateProvField Name of field in \code{df} with state/province names.
#' @param dfCountyField Name of field in \code{df} with county names.
#' @param dfSpeciesField Name of field in \code{df} with species names.
#' @examples
#' @export

matchDetectToGadm <- function(
	shape,
	df,
	species,
	dfStateProvField = 'stateProvince',
	dfCountyField = 'county',
	dfSpeciesField = 'species'
) {

	shape@data$detect <- 0

	if (!any(df[ , dfSpeciesField] != species)) {
	
		warning('Species not found in any records.')
		
	} else {
	
		df <- df[!is.na(df[ , dfSpeciesField]) & df[ , dfSpeciesField] == species, ]
	
		# lower case everything
		dfStateProv <- tolower(df[ , dfStateProvField])
		dfCounty <- tolower(df[ , dfCountyField])
		
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

		dfStateCounty <- paste(dfStateProv, dfCounty)
		shapeStateCounty <- paste(shapeStateProv, shapeCounty)

		df$FOUND <- FALSE
		
		# locate
		for (i in seq_along(shapeStateCounty)) {
		
			theseDfs <- which(dfStateCounty == shapeStateCounty[i])
			if (length(theseDfs) > 0) {
				shape@data$detect[i] <- length(theseDfs)
				df$FOUND[theseDfs] <- TRUE
			}
		
		}
		
		unfound <- sum(!df$FOUND)
		if (unfound > 0) warning(unfound, ' records were not matched in the shape.')
		
	}

	shape

}

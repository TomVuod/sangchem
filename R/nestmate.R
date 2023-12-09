#' @import dplyr
#' @export
find_sample_nestmate <- function(chr_ID, heterospecific=TRUE){
  if(is.na(chr_ID)) {
    warning("Chromatogram ID is NA")
    return(NA)
  }
  if(length(chr_ID)!=1) stop("'chr_ID' should be of length one.")
  data("development_data", envir = environment())
  if(!chr_ID %in% development_data$chromatogram_ID) stop(sprintf("Wrong chromatogram ID: %s", as.character(chr_ID)))
  focal_ant <- dplyr::filter(development_data, chromatogram_ID==chr_ID)
  partners <- dplyr::filter(development_data, colony==focal_ant$colony,
                            census_date==focal_ant$census_date, callow==0,
                            caste=="worker")
  if(heterospecific)
    partners <- dplyr::filter(partners, species!=focal_ant$species)
  else
    partners <- dplyr::filter(partners, species==focal_ant$species)
  if(nrow(partners)==0) return(NA)
  return(partners$chromatogram_ID)
}


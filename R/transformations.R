#' Peak proportions
#'
#' Represent MS data in the form of relative abundances of each peak
#'
#' @import dplyr
#' @export
relative_amounts <- function(chromatogram_data, chrID=NULL, peak_subset=NULL){
  if(!is.data.frame(chromatogram_data)) stop("'chromatogram_data' should be a data frame")
  if(nrow(chromatogram_data)==0) {
    warning("Empty data frame passed to 'relative_amounts'")
    return(NULL)
  }
  if (!missing(chrID) && any(!(chrID %in% chromatogram_data$chromatogram_ID))) {
    warning(sprintf("Chromatogram ID not found: %s", as.character(chrID)[!chrID %in% chromatogram_data$chromatogram_ID]))
    return(NULL)
  }
  if(missing(chrID)) {
    chrID <- unique(chromatogram_data$chromatogram_ID)
    if(length(chrID)!=1) stop("Guessing chromatogram ID from the data failed.")
  }
  dplyr::filter(chromatogram_data,chromatogram_ID %in% chrID,peak_ID>1) %>%
    dplyr::select(chromatogram_ID,peak_ID,abundance) -> filtered_data
  if (!is.null(peak_subset)) filtered_data <- dplyr::filter(filtered_data,peak_ID %in% peak_subset)
  if(nrow(filtered_data)==0) return(NULL)
  if(is.factor(filtered_data$chromatogram_ID)) filtered_data$chromatogram_ID <- droplevels(filtered_data$chromatogram_ID)
  filtered_data <- aggregate(list(abundance=filtered_data$abundance),
                             list(chromatogram_ID = filtered_data$chromatogram_ID, peak_ID=filtered_data$peak_ID), sum, drop=FALSE)
  filtered_data$abundance[is.na(filtered_data$abundance)] <- 0
  filtered_data <- split(filtered_data, filtered_data$chromatogram_ID) %>%
    lapply(function(x) {x$abundance <- x$abundance/sum(x$abundance); x}) %>%
    {function(x) purrr::reduce(x, merge, by = "peak_ID")}()
  filtered_data <- filtered_data[,!duplicated(colnames(filtered_data))] # remove extra peak ID columns
  filtered_data$relative_abundance <- apply(filtered_data[grep("abundance", colnames(filtered_data))], 1, function(x) mean(x))
  select(filtered_data, peak_ID, relative_abundance)
}

#' Relative abundances table
#'
#' Calculate relative abundances of peaks with chromatogram IDs as columns
#' @export
peak_proportions_table <- function(MS_data){
  if (!is.data.frame(MS_data)) stop("Data frame object missing")
  demanded_columns <- c('abundance','chromatogram_ID',"peak_ID")
  missing_columns <- setdiff(c('abundance','chromatogram_ID',"peak_ID"), colnames(MS_data))
  if (length(missing_columns)>0) stop(sprintf("Missing columns: %s", paste(missing_columns, collapse = ", ")))
  if(nrow(MS_data)==0) stop("Empty data frame")
  if (sum(is.na(MS_data[,"peak_ID"]))>0) warning(sprintf("Some rows contain NA values for the peak_ID"))
  MS_data <- MS_data[!is.na(MS_data[,"peak_ID"]),]
  chr_IDs <- c()
  prop_data <- lapply(split(MS_data, MS_data$chromatogram_ID), relative_amounts) %>%
    purrr::imap(function(x,y) {if(is.null(x)) return(NULL); colnames(x)[2] <- y; x}) %>%
    # account for NULL which might be returneb by 'relative_amounts'
    purrr::reduce(function(x,y,...) {if(is.null(y)) return(x); merge(x,y,...)}, all.x=TRUE, all.y=TRUE, by="peak_ID") %>%
    {function(x) {rownames(x) <- x[,1]; x[,-1]}}()
  prop_data[is.na(prop_data)] <- 0
  t(prop_data)
}

#' @export
clr_transformation <- function(x){
  log_geom_mean <- sum(log(x))/length(x)
  log(x) - log_geom_mean
}

#' Proportion of peak subset
#'
#' Calculate proportion of a set of peaks in samples
#'
#' @export
subsample_proportion<-function(peak_IDs, chromatogram_IDs, MS_data=NULL){
  if (is.null(peak_IDs)) stop("Peak IDs argument null")
  if (any(is.na(peak_IDs))) stop("A NA value in peak IDs argument")
  if (!all(is.numeric(peak_IDs))) stop("Peak IDs value should be numeric")
  if (is.null(chromatogram_IDs)) stop("Chromatogram IDs argument null")
  if (any(is.na(chromatogram_IDs))) stop("A NA value in Chromatogram IDs argument")
  if (!all(is.numeric(chromatogram_IDs))) warning("Chromatogram IDs are not numeric")
  if(missing(MS_data)){
    data("mass_spectra_data", envir=environment())
    MS_data <- mass_spectra_data
    rm(mass_spectra_data)
  }
  if (!all(chromatogram_IDs %in% unique(MS_data$chromatogram_ID)))
    warning("Not all chromatograms represented in the passed dataset")
  res<-c()
  for(i in 1:length(chromatogram_IDs)){
    prop_data <- relative_amounts(MS_data, chromatogram_IDs[i])
    res[i] <- sum(prop_data[prop_data$peak_ID %in% peak_IDs,"relative_abundance"])/sum(prop_data[,"relative_abundance"])
  }
  names(res)<-as.character(chromatogram_IDs)
  res
}

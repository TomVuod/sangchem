#' Peak proportions
#'
#' Represent MS data in the form of relative abundances of each peak
#'
#' @import dplyr
#' @export
relative_amounts<-function(chromatogram_data, chrID=NULL, peak_subset=NULL){
  if(!is.data.frame(chromatogram_data)) stop("'chromatogram_data' should be a data frame")
  if(nrow(chromatogram_data)==0) {
    warning("Empty data frame passed to 'relative_amounts'")
    return(NULL)
  }
  if (!missing(chrID) && !(is.numeric(chrID))) stop("Chromatogram ID should be numeric")
  if (!missing(chrID) && length(chrID)!=1) stop("One chromatogram ID should be passed in")
  if (!missing(chrID) && !(chrID %in% chromatogram_data$chromatogram_ID)) {
    warning(sprintf("Chromatogram ID not found: %s", as.character(chrID)))
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
  filtered_data <- aggregate(filtered_data$abundance,
                             list(peak_ID=filtered_data$peak_ID),sum)
  filtered_data <- rename(filtered_data,abundance=x)
  filtered_data <- filtered_data[,c("peak_ID", "abundance")]
  filtered_data$abundance<-filtered_data$abundance/sum(filtered_data$abundance)
  rename(filtered_data, relative_abundance=abundance)
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
  

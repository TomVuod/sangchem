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


#' Averaged CHC profile
#'
#' Construct CHC profile by averaging the input profiles
#'
#' @export
averaged_profile<-function(chr_IDs,MS_data,
                           CHC_amount_data, proportions=TRUE){
  if(length(chr_IDs)==0) stop("No chromatogram ID values in averaged_profile function")
  if(any(is.na(chr_IDs))) warning("NA chromatogram ID values in averaged_profile function")
  #if (length(chr_IDs)!=2) stop("Two chromatogram ID values should be passed in")
  output<-data.frame()
  for (i in 1:length(chr_IDs)){
    new_data<-relative_amounts(chromatogram_data= MS_data, chr_IDs[i])
    if (is.null(new_data)) next
    if (!proportions){
      mass=filter(CHC_amount_data,chromatogram_ID==chr_IDs[i])$corrected_mass
      if(is.na(mass)) return(NULL)
      new_data$relative_abundance=mass*new_data$relative_abundance
    }
    new_data$chromatogram_ID<-chr_IDs[i]
    output<-rbind(output,new_data)
  }
  if (nrow(output)==0) return(NULL)
  output<-aggregate(output$relative_abundance, list(output$peak_ID,
                                                    output$chromatogram_ID), sum, drop=FALSE)
  output$x[is.na(output$x)]<-0
  output<-rename(output,relative_abundance=x,peak_ID=Group.1,chromatogram_ID=Group.2)
  output<-aggregate(output$relative_abundance, list(output$peak_ID), mean)
  output=rename(output, relative_abundance=x, peak_ID=Group.1)
  if (!proportions){
    output=rename(output, abundance=relative_abundance)
  }
  output
}


average_by_sepcies <- function(dev_data, MS_data = MS_data, CHC_amounts = CHC_amounts){
  if (nrow(filter(dev_data,species=="F. sanguinea"))==0 | nrow(filter(dev_data,species=="F. fusca"))==0)
    return(NULL)
  average_fusca<-averaged_profile(filter(dev_data,species=="F. fusca")$chromatogram_ID,
                                  MS_data = MS_data,
                                  CHC_amount_data = CHC_amounts,
                                  proportions = FALSE)
  average_sanguinea<-averaged_profile(filter(dev_data,species=="F. sanguinea")$chromatogram_ID,
                                      MS_data = MS_data,
                                      CHC_amount_data = CHC_amounts,
                                      proportions = FALSE)
  if(is.null(average_fusca)|is.null(average_sanguinea)) return(NULL)
  joined<-full_join(average_fusca,average_sanguinea,by="peak_ID")
  joined$abundance.x[is.na(joined$abundance.x)]<-0
  joined$abundance.y[is.na(joined$abundance.y)]<-0
  joined
}

#' Average samples by sample
#'
#' Produce averaged samples by species and sampling date
#'
#' @export
sample_averaged<-function(dev_data, MS_data, CHC_amounts){
    dev_data %>%
    filter(callow==0,caste=="worker", chromatogram_ID %in% MS_data$chromatogram_ID) %>%
    group_by(colony, census_date) %>%
    group_map(~average_by_sepcies(.x, MS_data = mass_spectra_data, CHC_amounts = CHC_amounts)) %>%
    do.call(rbind, .)
}

#' Prportion of peak subset
#'
#' Calculate proportion of a set of peaks in samples
#'
#' @export
subsample_proportion<-function(peak_IDs, chromatogram_IDs, MS_data){
  if (is.null(peak_IDs)) stop("Peak IDs argument null")
  if (any(is.na(peak_IDs))) stop("A NA value in peak IDs argument")
  if (!all(is.numeric(peak_IDs))) stop("Peak IDs value should be numeric")
  if (is.null(chromatogram_IDs)) stop("Chromatogram IDs argument null")
  if (any(is.na(chromatogram_IDs))) stop("A NA value in Chromatogram IDs argument")
  if (!all(is.numeric(chromatogram_IDs))) warning("Chromatogram IDs are not numeric")
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

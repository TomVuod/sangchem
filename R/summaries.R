#' Peak proportions
#' 
#' Represent MS data in the form of relative abundances of each peak
#' 
#' @export
relative_amounts<-function(chrID, chromatogram_data, peak_subset=NULL){
  if (!(is.numeric(chrID))) stop("Chromatogram ID should be numeric")
  if (length(chrID)!=1) stop("One chromatogram ID should be passed in")
  if (!(chrID %in% chromatogram_data$chromatogram_ID)) {
    warning(sprintf("Chromatogram ID not found: %s", as.character(chrID)))
    return(NULL)
  }
  filter(chromatogram_data,chromatogram_ID %in% chrID,peak_ID>1) %>%
    select(chromatogram_ID,peak_ID,abundance) -> filtered_data
  if (!is.null(peak_subset)) filtered_data<-filter(filtered_data,peak_ID %in% peak_subset)
  filtered_data<-aggregate(filtered_data$abundance,
                           list(peak_ID=filtered_data$peak_ID),sum)
  filtered_data<-rename(filtered_data,abundance=x)
  filtered_data<-filtered_data[,c("abundance","peak_ID")]
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
    new_data<-relative_amounts(chr_IDs[i], chromatogram_data= MS_data)
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
    prop_data <- relative_amounts(chromatogram_IDs[i], MS_data)
    res[i] <- sum(prop_data[prop_data$peak_ID %in% peak_IDs,"relative_abundance"])/sum(prop_data[,"relative_abundance"])
  }
  names(res)<-as.character(chromatogram_IDs)
  res
} 
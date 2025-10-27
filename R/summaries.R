#' Averaged CHC profile
#'
#' Construct CHC profile by averaging the input profiles
#'
#' @export
averaged_profile<-function(chr_IDs,MS_data,
                           CHC_amount_data=NULL, proportions=TRUE){
  stopifnot(is.data.frame(MS_data))
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
  output <- aggregate(output$relative_abundance, list(output$peak_ID,
                                                      output$chromatogram_ID), sum, drop=FALSE)
  output$x[is.na(output$x)] <- 0
  output <- rename(output,relative_abundance=x,peak_ID=Group.1,chromatogram_ID=Group.2)
  mean_profile(output)
}

#' @export
mean_profile <- function(TIC_data, proportions=TRUE){
  output <- aggregate(TIC_data$relative_abundance, list(TIC_data$peak_ID), mean)
  output = rename(output, relative_abundance=x, peak_ID=Group.1)
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

#' Average samples
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

#' Calculate mean value for a selected variable per colony
#'
#' @export
colony_weighted_mean <- function(sample_data, var="normalized_mass"){
  if(missing(sample_data)){
    data("development_data", envir = environment())
    sample_data <- development_data
    rm(development_data)
  }
  sample_data %>% split(list(.$colony, .$census_date)) %>%
    lapply(function(x) data.frame(colony = x$colony[1], census_date=x$census_date[1],
                                  var=mean(x[[var]], na.rm=TRUE))) %>%
    {function(x) do.call(rbind, x)}() %>%
    split(.$colony) %>% lapply(function(x) data.frame(colony=x$colony[1], var = mean(x$var, na.rm=TRUE))) %>%
    {function(x) {df <- do.call(rbind, x); colnames(df)[2] <- var; df}}()
}

#' Calculate differential profile between callow and mature ants
#'
#' @export
differential_profile <- function(x){
  callow_IDs <- x[x$callow==1, "chromatogram_ID"]
  mature_IDs <- x[x$callow==0, "chromatogram_ID"]
  if(length(callow_IDs)==0 | length(mature_IDs)==0) return(NULL)
  average_profile_callow <- averaged_profile(callow_IDs, mass_spectra_data)
  average_profile_mature <- averaged_profile(mature_IDs, mass_spectra_data)
  diff_averaged_profile <- merge(average_profile_callow, average_profile_mature, all=TRUE, by="peak_ID")
  diff_averaged_profile[,"relative_abundance.x"][is.na(diff_averaged_profile[,"relative_abundance.x"])] <- 0
  diff_averaged_profile[,"relative_abundance.y"][is.na(diff_averaged_profile[,"relative_abundance.y"])] <- 0
  data.frame(peak_ID = diff_averaged_profile[,"peak_ID"], relative_abundance = diff_averaged_profile[,"relative_abundance.x"] - diff_averaged_profile[,"relative_abundance.y"])
}

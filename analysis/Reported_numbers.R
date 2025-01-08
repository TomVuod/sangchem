# This script caontains the way of the calculation of numbers reported in the manuscript

library(devtools)
# Total number of peaks
data("mass_spectra_data")
data("development_data")
data("separation_experiment_data")
length(setdiff(unique(mass_spectra_data[mass_spectra_data$chromatogram_ID %in% c(development_data$chromatogram_ID,
                                                                         separation_experiment_data$chromatogram_ID), "peak_ID"]), 0))

# number of F. sanguinea markers
peak_sang <- intersect(sangchem:::peak_sang, mass_spectra_data$peak_ID)
length(peak_sang)

# mean mass of the F. sanguinea markers
data("field_colonies")
library(stringi)
sanguinea_samples <- as.character(unique(field_colonies$chromatogram_ID)[stri_detect_regex(unique(field_colonies$chromatogram_ID),"sangFs[0-9].?\\([0-9]\\)")])
fusca_samples <- as.character(unique(field_colonies$chromatogram_ID)[stri_detect_regex(unique(field_colonies$chromatogram_ID), "fus[0-9].?-[0-9]\\([1-2]\\)")])

subsample_proportion <- function(peak_IDs, chromatogram_IDs, MS_data=NULL){
  if (is.null(peak_IDs)) stop("Peak IDs argument null")
  if (any(is.na(peak_IDs))) stop("A NA value in peak IDs argument")
  if (!all(is.numeric(peak_IDs))) stop("Peak IDs value should be numeric")
  if (is.null(chromatogram_IDs)) stop("Chromatogram IDs argument null")
  if (any(is.na(chromatogram_IDs))) stop("A NA value in Chromatogram IDs argument")
  if(missing(MS_data)){
    data("mass_spectra_data", envir=environment())
    MS_data <- mass_spectra_data
    rm(mass_spectra_data)
  }
  if (!all(chromatogram_IDs %in% unique(MS_data$chromatogram_ID)))
    warning("Not all chromatograms represented in the passed dataset")
  res<-c()
  for(i in 1:length(chromatogram_IDs)){
    prop_data <- sangchem:::relative_amounts(MS_data, chromatogram_IDs[i])
    res[i] <- sum(prop_data[prop_data$peak_ID %in% peak_IDs,"relative_abundance"])/sum(prop_data[,"relative_abundance"])
  }
  names(res)<-as.character(chromatogram_IDs)
  res
}

mean(subsample_proportion(peak_sang, sanguinea_samples, field_colonies)) # in F. sanguinea
mean(subsample_proportion(peak_sang, fusca_samples, field_colonies)) # in F. fusca


# number of F. fusca markers
peak_fusca <- intersect(sangchem:::peak_fusca, mass_spectra_data$peak_ID)
length(peak_fusca)

# mean mass of the F. fusca markers

mean(subsample_proportion(peak_fusca, sanguinea_samples, field_colonies)) # in F. sanguinea
mean(subsample_proportion(peak_fusca, fusca_samples, field_colonies)) # in F. fusca

# number of callow markers
length(sangchem:::peaks_callow)


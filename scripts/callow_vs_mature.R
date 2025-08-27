compare_age_groups <- function(x){
  x$peak_number <- unlist(lapply(x$chromatogram_ID, function(x) nrow(mass_spectra_data[mass_spectra_data$chromatogram_ID==x & mass_spectra_data$peak_ID > 1,])))
  x
}

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


suppressMessages({
  library(sangchem)
  library(dplyr)
  library(mixOmics)
})

# load data
data("development_data")
data("mass_spectra_data")

DA_data <- dplyr::filter(development_data, caste=="worker", !remarks %in% c("aggression actor", "aggression target", "queenless colony"), species=="F. sanguinea") %>%
  dplyr::select(colony, species, callow, census_date, chromatogram_ID, sang_prop) %>%
  split(list(.$colony, .$census_date), drop=TRUE) %>%
  lapply(function(x){if(sum(x$callow)==0) return(NULL); x}) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  lapply(differential_profile) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  purrr::reduce(function(x,y) merge(x,y,all=TRUE,by="peak_ID"))

DA_data <- as.data.frame(lapply(as.list(DA_data), function(x) {x[is.na(x)] <- 0; x}))

df <- data.frame()
for(i in 1:nrow(DA_data)) df <- rbind(df, data.frame(peak_ID = DA_data[i,"peak_ID"], mean_change = mean(as.numeric(DA_data[i,2:ncol(DA_data)]))))


boxplot(DA_data$peak_number~DA_data$callow)


age_profiles <- dplyr::filter(development_data, caste=="worker", !remarks %in% c("aggression actor", "aggression target", "queenless colony"), species=="F. sanguinea") %>%
  dplyr::select(colony, species, callow, census_date, chromatogram_ID, sang_prop) %>%
  split(list(.$colony, .$census_date), drop=TRUE) %>%
  lapply(function(x){if(sum(x$callow)==0) return(NULL); x}) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  lapply(age_profile_sample) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  split(.$callow, drop=TRUE) %>%
  lapply(combine_profile)


age_profile_sample <- function(x){
  callow_IDs <- x[x$callow==1, "chromatogram_ID"]
  mature_IDs <- x[x$callow==0, "chromatogram_ID"]
  if(length(callow_IDs)==0 | length(mature_IDs)==0) return(NULL)
  average_profile_callow <- averaged_profile(callow_IDs, mass_spectra_data)
  average_profile_mature <- averaged_profile(mature_IDs, mass_spectra_data)
  average_profile_callow$callow <- 1
  average_profile_mature$callow <- 0
  rbind(average_profile_callow,  average_profile_mature)
}




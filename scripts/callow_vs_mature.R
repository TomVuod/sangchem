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

colnames(DA_data)[2:length(colnames(DA_data))] <- paste0("col_", 1:(length(colnames(DA_data))-1))
diff_data <- tidyr::gather(DA_data, key = "sample", value = "relative_abundance", -peak_ID)

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



age_profiles <- dplyr::filter(development_data, caste=="worker", !remarks %in% c("aggression actor", "aggression target", "queenless colony"), species=="F. sanguinea") %>%
  dplyr::select(colony, species, callow, census_date, chromatogram_ID, sang_prop) %>%
  split(list(.$colony, .$census_date), drop=TRUE) %>%
  lapply(function(x){if(sum(x$callow)==0) return(NULL); x}) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  lapply(age_profile_sample) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  purrr::reduce(rbind) %>%
  split(.$callow, drop=TRUE) %>%
  lapply(mean_profile)


all_peaks <- unique(c(age_profiles[[1]]$peak_ID, age_profiles[[2]]$peak_ID))
missing_peaks <- setdiff(all_peaks, age_profiles[[1]]$peak_ID)
age_profiles[[1]] <- rbind(age_profiles[[1]], data.frame(peak_ID = missing_peaks, relative_abundance = 0))
age_profiles[[1]]$age = "mature"
missing_peaks <- setdiff(all_peaks, age_profiles[[2]]$peak_ID)
age_profiles[[2]] <- rbind(age_profiles[[2]], data.frame(peak_ID = missing_peaks, relative_abundance = 0))
age_profiles[[2]]$age = "callow"
#plot_data <- rbind(age_profiles[[1]], age_profiles[[2]])

library(ggplot2)

ggplot(NULL, aes(x=peak_ID, y=relative_abundance, fill = age, color=age))+
  geom_bar(data = age_profiles[[1]], stat = "identity",linewidth=0.2,
           position = position_dodge(0.9), width = 0.4)+
  geom_bar(data = age_profiles[[2]], stat = "identity",linewidth=0.2,
           position = position_dodge(0.9), width = 0.4)+
  scale_color_manual(values=c("transparent", "#00000099"))+
  scale_fill_manual(values=c("#99999966", "transparent"))+
  geom_text(data = age_profiles[[1]], aes(label=peak_ID), size = 2,
             hjust=-0.8, color="black", angle = 90, vjust = 0.3)+
  geom_point(data = diff_data, aes(x=peak_ID, y=relative_abundance), color = "black", fill="black", size = 0.1)+
  theme(panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_rect(color="white", fill = "white"))



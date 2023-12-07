# change to development_data

development_data[development_data$colony=="SD18-11"&development_data$sang_prop>0.1,"remarks"] <- "queenless colony"

separation_experiment_data <- separation_experiment_data[,!colnames(separation_experiment_data) %in%
                                                           c("X", "census_date", "pupae_number", "species")]

separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"end"] <-
  as.POSIXct(separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"end"]) + lubridate::years(1)

separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"delta_1"] <-
  as.POSIXct(separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"end"]) -
  separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"pupa_1"]

separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"delta_2"] <-
  as.POSIXct(separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"end"]) -
  separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"pupa_2"]


separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"mean_delta"] <-
  mean(separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"delta_1"],
       separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"delta_2"])

# remove sample with poor chromarogram
separation_experiment_data <- separation_experiment_data[separation_experiment_data$chromatogram_ID!=442,]

mass_spectra_data[mass_spectra_data$peak_ID==7, "identification"] <- "C23"
mass_spectra_data[mass_spectra_data$peak_ID==14, "identification"] <- "C24"
mass_spectra_data[mass_spectra_data$peak_ID==24, "identification"] <- "C25"
mass_spectra_data[mass_spectra_data$peak_ID==35, "identification"] <- "C26"
mass_spectra_data[mass_spectra_data$peak_ID==47, "identification"] <- "C27"
mass_spectra_data[mass_spectra_data$peak_ID==56, "identification"] <- "C28"
mass_spectra_data[mass_spectra_data$peak_ID==68, "identification"] <- "C29"
mass_spectra_data[mass_spectra_data$peak_ID==81, "identification"] <- "C30"
mass_spectra_data[mass_spectra_data$peak_ID==90, "identification"] <- "C31"


separation_experiment_data$colony[!separation_experiment_data$colony %in% c("17-1", "17-2")] <-
  paste0("SD",separation_experiment_data$colony[!separation_experiment_data$colony %in% c("17-1", "17-2")])

separation_experiment_data$colony[separation_experiment_data$colony %in% c("17-1", "17-2")] <-
  paste0("W",separation_experiment_data$colony[separation_experiment_data$colony %in% c("17-1", "17-2")])


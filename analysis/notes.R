# change to development_data

development_data[development_data$colony=="SD18-11"&development_data$sang_prop>0.1,"remarks"] <- "queenless colony"

separation_experiment_data <- separation_experiment_data[,!colnames(separation_experiment_data) %in%
                                                           c("X", "census_date", "pupae_number", "species")]

separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"end"] <-
  as.POSIXct(separation_experiment_data[separation_experiment_data$chromatogram_ID == 446,"end"]) + years(1)

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

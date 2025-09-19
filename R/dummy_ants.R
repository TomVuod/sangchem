recognize_groups<-function(sample_data=separation_experiment_data){
  filter(sample_data,treatment_id>4) %>% group_by(colony, year) -> temp
  groups<-group_indices(temp)
  temp<-ungroup(temp)
  temp$group<-groups
  temp
}

determine_donor_chr_ID<-function(chrID, species, donor_mapping,
                                 sep_data){
  col <- as.character(filter(sep_data,chromatogram_ID==chrID)$colony)
  year <- filter(sep_data,chromatogram_ID==chrID)$year
  if((species=="F. fusca")&(col=="SD18-5")&(year==2020)){
    pupa_1<-filter(sep_data,chromatogram_ID==chrID)$pupa_1
    if(pupa_1<as.POSIXct("2020-08-01",format="%Y-%m-%d"))
      return(293)
    else return(294)
  }
  species_ = species
  year_ = year
  donor_chr_ID <- filter(donor_mapping,species==species_,colony==col,year==year_)$chromatogram_ID
  if(length(donor_chr_ID)==0) return(NA)
  else donor_chr_ID
}

# add control sample if present; if not removes the treatment sample
pair_present<-function(sample_data,treatment_id_){
  if (sum(sample_data$treatment_id==treatment_id_)==0)
    return(NULL)
  treatments<-5:8
  control_id<-rep(c(9,10),2)[treatments==treatment_id_]
  if (sum(sample_data$treatment_id==control_id)==0)
    return(NULL)
  filter(sample_data,treatment_id %in% c(treatment_id_,control_id))
}


select_group_with_control<-function(sample_data,treatment_id){
  sample_data<-split(sample_data,sample_data$group)
  sample_data<-lapply(sample_data, function(x) pair_present(x,treatment_id_=treatment_id))
  do.call(rbind,sample_data)
}

# calculate distance to the colony which was used as donor of CHCs applicated onto the dummy ants
# the same is done for corresponding control samples
calculate_distance_to_donor <- function(treatment_id, sample_data, MS_data, donor_mapping){
  sample_data <- recognize_groups(sample_data)
  sample_data <- select_group_with_control(sample_data,treatment_id)
  treatments <- 5:8
  species<-rep(c("F. fusca", "F. sanguinea"),each=2)[treatments==treatment_id]
  sample_data$donor_id<-mapply(function(x) determine_donor_chr_ID(x, species,
                                                                  sep_data = sample_data, donor_mapping=donor_mapping), sample_data$chromatogram_ID)
  sample_data <- filter(sample_data, !is.na(donor_id))
  sample_data$distance<-mapply(function(x,y) pair_distance(x,y, MS_data), sample_data$chromatogram_ID, sample_data$donor_id)
  sample_data %>% group_by(colony,year,group,treatment) %>% summarise(mean_dist=mean(distance))
}

calculate_C22_fraction<-function(treatment_id, sample_data, MS_data){
  sample_data <- recognize_groups(sample_data)
  sample_data<-recognize_groups(sample_data)
  sample_data<-select_group_with_control(sample_data,treatment_id)
  sample_data$C22_prop<-mapply(function(x, y, peak_areas) peak_proportion(y, x, peak_areas),
                               sample_data$chromatogram_ID, MoreArgs = list(y=1, peak_areas=MS_data))
  sample_data %>% group_by(colony,year,group,treatment) %>% summarise(C22_prop=mean(C22_prop))
}

#' Dummny ant effect
#'
#' Test of the effect of dummy ant on the CHC profile of separated workers
#'
#' @param treatment_id a numeric indicating treatment ID to be tested on the effect of
#' dummy ant.
#' @param sample_data a data frame describing samples used in separation experiment; retrieved with `data("separation_experiment_data")`
#' @param MS_data a data frame with the results of the chromatogram peaks integration; ; retrieved with `data("mass_spectra_data")`
#' @param donor_mapping a data frame providing sample ID mapping between control and treatment samples
#' in the dummy ant experiment
#' @export
test_dummy_ant_effect <- function(treatment_id, sample_data, MS_data, donor_mapping,
                                  caption=""){
  distances_to_donor<-calculate_distance_to_donor(treatment_id = treatment_id, sample_data = sample_data,
                                                  donor_mapping=donor_mapping, MS_data=MS_data)
  split_by_treatment <- split(distances_to_donor, distances_to_donor$treatment)
  to_test <- full_join(split_by_treatment[[2]], split_by_treatment[[1]], by=c("group","colony","year"))
  to_test <- to_test[,c("colony", "treatment.x", "mean_dist.x", "treatment.y", "mean_dist.y")]
  to_test[,c("mean_dist.x", "mean_dist.y")] <- round(to_test[,c("mean_dist.x", "mean_dist.y")], 3)
  args <- list(x=to_test, caption = caption,
               col.names = c("Colony ID", "Treamtent", "Averaged distance", "Contrast", "Averaged distance"))
  print(wilcox.test(to_test$mean_dist.x,to_test$mean_dist.y,paired = TRUE))
  return(args)
}


#' @export
test_C22_proportion <- function(treatment_id, sample_data, MS_data,
                                caption=""){
  prop_data <- calculate_C22_fraction(treatment_id, sample_data, MS_data)
  split_by_treatment <- split(prop_data, prop_data$treatment)
  to_test<-full_join(split_by_treatment[[2]], split_by_treatment[[1]], by=c("group","colony","year"))
  to_test <- to_test[,c("colony", "treatment.x", "C22_prop.x", "treatment.y", "C22_prop.y")]
  to_test[,c("C22_prop.x", "C22_prop.y")] <- round(to_test[,c("C22_prop.x", "C22_prop.y")], 3)
  args <- list(x=to_test, caption = caption,
               col.names = c("Colony ID", "Treamtent", "Averaged C22 fraction", "Contrast", "Averaged C22 fraction"))
  print(wilcox.test(to_test$C22_prop.x,to_test$C22_prop.y,paired = TRUE))
  return(args)
}

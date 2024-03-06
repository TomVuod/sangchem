recognize_groups<-function(sample_data=separation_experiment_data){
  filter(sample_data,treatment_id>4) %>% group_by(colony, year) -> temp
  groups<-group_indices(temp)
  temp<-ungroup(temp)
  temp$group<-groups
  temp
}

determine_donor_chr_ID<-function(chrID,species_,donor_mapping_=donor_mapping,
                                 sep_data=separation_experiment_data){
  col<-as.character(filter(sep_data,chromatogram_ID==chrID)$colony)
  year_<-filter(sep_data,chromatogram_ID==chrID)$year
  if((species_=="F. fusca")&(col=="SD18-5")&(year_==2020)){
    pupa_1<-filter(sep_data,chromatogram_ID==chrID)$pupa_1
    if(pupa_1<as.POSIXct("2020-08-01",format="%Y-%m-%d"))
      return(293)
    else return(294)
  }
  filter(donor_mapping_,species==species_,colony==col,year==year_)$chromatogram_ID
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
calculate_distance_to_donor<-function(treatment_id,sample_data){
  sample_data<-recognize_groups(sample_data)
  sample_data<-select_group_with_control(sample_data,treatment_id)
  treatments<-5:8
  species<-rep(c("F. fusca", "F. sanguinea"),each=2)[treatments==treatment_id]
  sample_data$donor_id<-mapply(function(x) determine_donor_chr_ID(x,species), data$chromatogram_ID)
  sample_data$distance<-mapply(function(x,y) pair_distance(c(x,y)), sample_data$chromatogram_ID, sample_data$donor_id)
  sample_data %>% group_by(colony,year,group,treatment_id) %>% summarise(mean_dist=mean(distance))
}

#' Dummny ant effect
#'
#' Test of the effect of dummy ant on the CHC profile of separated workers
#'
#' @param treatment_id a numeric indicating treatment ID to be tested on the effect of
#' @param sample_data a data frame describing samples used in separation experiment; rerieved with `data("separation_experiment_data")`
#' dummy ant.
#' @export
test_dummy_ant_effect <- function(treatment_id, sample_data){
  distances_to_donor<-calculate_distance_to_donor(treatment_id = treatment_id, sample_data = sample_data)
  splitted<-split(distances_to_donor, distances_to_donor$treatment_id)
  to_test<-full_join(splitted[[1]], splitted[[2]], by=c("group","colony","year"))
  print(to_test)
  print(wilcox.test(to_test$mean_dist.x,to_test$mean_dist.y,paired = TRUE))
}



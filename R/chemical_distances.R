global_distances_matrix <- function(){
  data("separation_experiment_data", envir = environment())
  data("development_data", envir = environment())
  distances_matrix <- matrix(NA, nrow = length(unique(separation_experiment_data$chromatogram_ID)),
                             ncol = length(unique(na.omit(colony_to_chrID))),
                             dimnames = list(unique(separation_experiment_data$chromatogram_ID),
                                             unique(na.omit(colony_to_chrID))))
  set <- function(sepID, devID, val) {
    distances_matrix[as.character(sepID), as.character(devID)] <<- val
  }
  get <- function(sepID, devID) {
    distances_matrix[as.character(sepID), as.character(devID)]
  }
  return(list(set=set, get=get))
}

#' @importFrom vegan vegdist
#' @export
pair_distance <- function(chr_ID_1, chr_ID_2, MS_data, distances_registry=NULL,
                          ...){
  if (any(is.na(c(chr_ID_1, chr_ID_2)))) return(NA)
  stopifnot((length(chr_ID_1)==1) & (length(chr_ID_2)==1))
  if (!is.null(distances_registry) && !is.na(distances_registry$get(chr_ID_1, chr_ID_2)))
    return(distances_registry$get(chr_ID_1, chr_ID_2))
  data_table <- full_join(relative_amounts(MS_data[MS_data$chromatogram_ID==chr_ID_1,]),
                          relative_amounts(MS_data[MS_data$chromatogram_ID==chr_ID_2,]),
                        by="peak_ID")
  data_table<-as.matrix(select(data_table,-peak_ID))
  data_table[is.na(data_table)] <- 0
  dist <- vegdist(t(data_table), method="bray")[[1]]
  if(!is.null(distances_registry)) distances_registry$set(chr_ID_1, chr_ID_2, dist)
  dist
}

distances_to_source_colonies <- function(GCMS_samples, MS_data, colony_to_chrID,
                                      slave_source_colony, randomize=FALSE,
                                      distances_registry=NULL){
  GCMS_samples$chr_ID_1 <- GCMS_samples$chromatogram_ID
  if(randomize){
    GCMS_samples$colony <- factor(GCMS_samples$colony)
    levels(GCMS_samples$colony) <- sample(levels(GCMS_samples$colony))
    GCMS_samples$colony <- as.character(GCMS_samples$colony)
  }
  GCMS_samples$chr_ID_2 <- colony_to_chrID[slave_source_colony[GCMS_samples$colony]]
  GCMS_samples <- GCMS_samples[!is.na(GCMS_samples$chr_ID_2),]
  GCMS_samples$dissimilarity <- purrr::pmap_dbl(GCMS_samples, pair_distance, MS_data=MS_data,
                                                distances_registry=distances_registry)
  mean(colony_mean_distance(GCMS_samples)$dissimilarity)
}

empirical_p_value<-function(val,distribution){
  if(any(is.na(val))) warning("NA values in the treatment data")
  if(length(val)==0) stop("No treatment data")
  if(length(val)==1)
    sum(na.omit(distribution)<=val)/length(na.omit(distribution))
}

#' @export
distances_permutation_test <- function(treatment_id, n_iter=1e5+1, exclude_naked_pupae=FALSE){
  distances_registry <- global_distances_matrix()
  data("separation_experiment_data", envir=environment())
  data("mass_spectra_data", envir=environment())
  data("colony_to_chrID", envir=environment())
  data("slave_source_colony", envir=environment())
  separation_data <- separation_experiment_data[separation_experiment_data$treatment_id == treatment_id,]
  if(exclude_naked_pupae) separation_data <- separation_data[apply(separation_data[c("cocoon_1", "cocoon_2")], 1, function(x) sum(as.numeric(x))==2),]
  true_distance <- distances_to_source_colonies(separation_data, mass_spectra_data, colony_to_chrID,
                                                 slave_source_colony, distances_registry = distances_registry)
  randomized_distances <- unlist(lapply(seq_len(n_iter), function(x) distances_to_source_colonies(separation_data,
                                                                                       mass_spectra_data,
                                                                                       colony_to_chrID,
                                                                                       slave_source_colony,
                                                                                       randomize = TRUE,
                                                                                       distances_registry = distances_registry)))
 list(p.val = empirical_p_value(true_distance, randomized_distances),
      true_distance=true_distance, rand_distances = randomized_distances)
}





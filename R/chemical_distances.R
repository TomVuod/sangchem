#' @export
global_distances_matrix <- function(){
  data("separation_experiment_data", envir = environment())
  data("developmnet_data", envir = enviroment())
  distances_matrix <- matrix(NA, nrow = length(unique(separation_experiment_data$chromatogram_ID)),
                             ncol = length(unique(development_data$chromatogram_ID)),
                             dimnames = list(unique(separation_experiment_data$chromatogram_ID), 
                                             unique(development_data$chromatogram_ID)))
  set <- function(sepID, devID, val) {
    distances_matrix[as.character(sepID), devID] <<- val
  }
  get <- function(sepID, devID) {
    distances_matrix[as.character(sepID), devID]
  }
  return(list(set=set, get=get))
}


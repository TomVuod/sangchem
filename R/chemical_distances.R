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
  clear <- function()
    distances_matrix[!is.na(distances_matrix)] <<- NA
}

#' @importFrom vegan vegdist
#' @export
pair_distance <- function(chr_ID_sep, chr_ID_source, MS_data){
  if (any(is.na(c(chr_ID_sep, chr_ID_source)))) return(NA)
  data_table <- full_join(relative_amounts(chr_ID_sep),relative_amounts(chr_ID_source),
                        by="peak_ID")
  data_table<-as.matrix(select(data_table,-peak_ID))
  data_table[is.na(data_table)] <- 0
  vegdist(t(data_table),method="bray")[[1]]
}

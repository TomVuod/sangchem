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

<<<<<<< HEAD
#' Vartiable selection based on their quality
#'
#' Select peaks based on their frequency across the samples
#' and mean abundance (from non-zero entries)
#'
#' @export
freq_abun_QC <- function(prop_data, freq_cutoff = 0.1, abundance_cutoff = 0.002){
  fr_select <- apply(prop_data, 2, function(x) sum(x>0)/length(x)) >= freq_cutoff
  prop_data <- prop_data[,fr_select]
  # re-normalize to 0
  prop_data <- prop_data/rowSums(prop_data)
  mean_abundances <- apply(prop_data, 2, function(x) mean(x[x>0]))
  while(any(mean_abundances < abundance_cutoff)){
    lowest_peak = which.min(mean_abundances)
    prop_data <- prop_data[,-lowest_peak]
    prop_data <- prop_data/rowSums(prop_data)
    mean_abundances <- apply(prop_data, 2, function(x) mean(x[x>0]))
  }
  prop_data
}
=======
#' Vartiable selection based on their quality
#' 
#' Select peaks based on their frequency across the samples
#' and mean abundance (from non-zero entries)
#' 
#' @export
freq_abun_QC <- function(prop_data, freq_cutoff = 0.1, abundance_cutoff = 0.002){
  fr_select <- apply(prop_data, 2, function(x) sum(x>0)/length(x)) >= freq_cutoff
  abund_select <- apply(prop_data, 2, function(x) mean(x[x>0])) >= abundance_cutoff
  prop_data[,fr_select&abund_select]
}
>>>>>>> supp

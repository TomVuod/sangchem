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
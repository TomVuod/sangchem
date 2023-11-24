#' Calculate species identity index
#'
#' @export
calculate_SII <- function(MS_data, prediction_model=species_prediction_model,
                          PCA_res = PCA_species_discr){
  peak_prop <- peak_proportions_table(MS_data)
  peak_prop <- peak_prop[,colnames(peak_prop) %in% rownames(PCA_res$rotation)]
  peak_prop <- peak_prop/rowSums(peak_prop)
  missing_variables <- setdiff(rownames(PCA_res$rotation), colnames(peak_prop))
  peak_prop <- cbind(peak_prop, matrix(0, nrow=nrow(peak_prop), ncol=length(missing_variables),
                                       dimnames = list(rownames(peak_prop), missing_variables)))
  peak_prop <- peak_prop[,rownames(PCA_res$rotation)]

  # PCA transformation using data from training set
  peak_prop <- t(apply(peak_prop, 1, function(x) x-PCA_res$center))
  peak_prop <- peak_prop%*%PCA_res$rotation
  predicted <- predict(prediction_model, peak_prop)$predict
  cbind(data.frame(chromatogram_ID=as.numeric(dimnames(predicted)[[1]])),
        predicted_species=predicted[,2,prediction_model$ncomp])
}

#' @export
find_nestmates_SII <- function(chr_IDs, species_indices, heterospecific=TRUE){
  lapply(chr_IDs, function(x) find_sample_nestmate(x, heterospecifc=heterospecific)) %>%
  lapply(function(x) species_indices[species_indices$chromatogram_ID %in% x,"predicted_species"]) %>%
  lapply(function(x) mean(x, na.rm=TRUE)) %>%
  unlist()
}

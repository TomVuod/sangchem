#' Calculate species identity index
#'
#' @export
calculate_SII <- function(normalized_TIC, PLS_DA_model=species_prediction_model,
                          PCA_res = PCA_species_discr, predicted_category=2,
                          multilevel=NULL){
  normalized_TIC <- normalized_TIC[,colnames(normalized_TIC) %in% rownames(PCA_res$rotation),drop=FALSE]
  normalized_TIC <- normalized_TIC/rowSums(normalized_TIC)
  missing_variables <- setdiff(rownames(PCA_res$rotation), colnames(normalized_TIC))
  normalized_TIC <- cbind(normalized_TIC, matrix(0, nrow=nrow(normalized_TIC), ncol=length(missing_variables),
                                                   dimnames = list(rownames(normalized_TIC), missing_variables)))
  normalized_TIC <- normalized_TIC[,rownames(PCA_res$rotation),drop=FALSE]
  normalized_TIC[normalized_TIC==0] <- 10^-16
  transformed_TIC <- t(apply(normalized_TIC,1,clr_transformation))
  # PCA transformation using data from training set
  transformed_TIC_PC <- predict(PCA_res, transformed_TIC)
  if(!is.null(multilevel)) transformed_TIC_PC<-withinVariation(transformed_TIC_PC, data.frame(multilevel))
  predicted <- predict(PLS_DA_model, transformed_TIC_PC)$predict
  cbind(data.frame(chromatogram_ID=as.numeric(dimnames(predicted)[[1]])),
        predicted_species=predicted[,predicted_category,PLS_DA_model$ncomp])
}

#' @export
find_nestmates_SII <- function(chr_IDs, species_indices, heterospecific=TRUE){
  lapply(chr_IDs, function(x) find_sample_nestmate(x, heterospecific=heterospecific)) %>%
  lapply(function(x) species_indices[species_indices$chromatogram_ID %in% x,"predicted_species"]) %>%
  lapply(function(x) mean(x, na.rm=TRUE)) %>%
  unlist()
}

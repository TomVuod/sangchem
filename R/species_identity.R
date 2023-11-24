#' Calculate species identity index
#'
#' @export
calculate_SII <- function(MS_data, prediction_model=species_prediction_model,
                          PCA_res = PCA_species_discr){
  peak_prop <- peak_proportions_table(MS_data)
  peak_prop <- peak_prop[,peak_prop %in% rownames(PCA_res$rotation)]
  peak_prop <- peak_prop/rowSums(peak_prop)
  missing_variables <- setdiff(rownames(PCA_res$rotation), colnames(peak_prop))
  peak_prop <- cbind(peak_prop, matrix(0, nrow=nrow(peak_prop), ncol=length(missing_variables),
                                       dimnames = list(rownames(peak_prop), missing_variables)))
  peak_prop <- peak_prop[,rownames(PCA_res$rotation)]

  # PCA transformation using data from training set
  peak_prop <- t(apply(peak_prop, 1, function(x) x-PCA_res$center))
  peak_prop <- peak_prop%*%PCA_res$rotation
  predicted <- predict(species_model, peak_prop)$predict
  cbind(data.frame(chromatogram_ID=as.numeric(dimnames(predicted)[[1]])),
        predicted_species=predicted[,1,species_model$ncomp])
}




predict_species<-function(chrIDs, model=plsda_CBII, ncomp_included=3){
  if (any(is.na(chrIDs))) {
    warning("NA value of chromatogram ID passed in")
    return(NA)
  }
  data<-pca_projection(chrIDs)
  predicted<-predict(model,data)$predict
  cbind(data.frame(chromatogram_ID=as.numeric(dimnames(predicted)[[1]])),
        predicted_species=predicted[,1,ncomp_included])
}

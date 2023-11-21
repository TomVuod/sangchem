run_linear_model <-function (averaged_data){
  if (!is.data.frame(averaged_data)) stop("'averaged_data' should be a data frame.")
  if (!"abundance.x" %in% colnames(averaged_data))  stop("'abundance.x' column absent.")
  if (!"abundance.y" %in% colnames(averaged_data))  stop("'abundance.y' column absent.")
  if (nrow(averaged_data)==0) {
    warning("Empty data frame passed to 'run_linear_model function")
    return(NULL)
  }
  model_res<-lm(abundance.x~abundance.y,data=averaged_data)
  n=nrow(averaged_data)
  correlation=cor(averaged_data$abundance.x, averaged_data$abundance.y)
  p_val <- ifelse(nrow(model_res$model)>1,summary(model_res)$coeff[2,4], NA)
  data.frame(peak_ID=averaged_data$peak_ID[1], correlation=correlation, p_val=p_val, n=n)
}

#' Correlations in CHC change
#'
#' Calculate how CHC amounts across samples correlate between heterospecific nestamates
#'
#' @export
peak_correlations <- function(averaged_data){
  averaged_data <- split(averaged_data, averaged_data$peak_ID)
  averaged_data <- lapply(averaged_data, run_linear_model)
  do.call(rbind, averaged_data)
}

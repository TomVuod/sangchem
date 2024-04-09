library(sangchem)
library(mixOmics)
library(dplyr)

calculate_AUC <- function(x,y){
  y <- y[order(x)]
  x <- x[order(x)]
  non_unique_x_ind <- which(duplicated(x))
  non_unique_x_ind <- sort(unique(c(non_unique_x_ind, non_unique_x_ind-1)))
  y[non_unique_x_ind] <- sort(y[non_unique_x_ind])
  if(any(diff(y)<0)) stop("y should be non-decreasing function of x")
  x_intervals <- diff(x)
  pair_mean_y <- (y[1:(length(y)-1)] + y[2:length(y)])/2
  sum(x_intervals*pair_mean_y)
}


calculate_accuracy_metrics <- function(values, labels){
  labels <- labels[order(values)]
  values <- sort(values)
  res_list <- list()
  all_combinations <- data.frame(threshold_reached=as.logical(c(1,1,0,0)), category = as.logical(c(1,0,1,0)))
  confusion_matrix <- data.frame(threshold_reached =0, category = as.numeric(labels))
  confusion_matrix <- rbind(all_combinations, confusion_matrix)
  tab <- table(confusion_matrix)
  tab <- tab - 1
  res_list[[1]] <- c(TP = tab[2,2], FP = tab[2, 1], TN = tab[1, 1],
                     FN = tab[1, 2], cutoff = NA)
  for(i in seq_along(values)){
    threshold_reached <- values >= values[i]
    confusion_matrix <- data.frame(threshold_reached =threshold_reached, category = as.numeric(labels))
    confusion_matrix <- rbind(all_combinations, confusion_matrix)
    tab <- table(confusion_matrix)
    # account for adding one observation for combination
    tab <- tab - 1
    res_list[[i+1]] <- c(TP = tab[2,2], FP = tab[2, 1], TN = tab[1, 1],
                         FN = tab[1, 2], cutoff = values[i])
  }
  res <- as.data.frame(do.call(rbind, res_list))
  TPR <- res$TP/(res$TP + res$FN)
  FPR <- res$FP/(res$FP + res$TN)
  list(TPR = TPR, FPR = FPR, cutoff = res$cutoff, confusion_matrix_data = res)
}
# load data
data("development_data")
data("mass_spectra_data")

DA_data <- dplyr::filter(development_data, caste=="worker", !remarks %in% c("aggression actor", "aggression target", "queenless colony")) %>%
  dplyr::select(colony, species, callow, census_date, chromatogram_ID, sang_prop) %>%
  split(list(.$colony, .$census_date), drop=TRUE) %>%
  lapply(function(x){if(sum(x$callow)==0) return(NULL); x}) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  purrr::map2(seq_along(.), function(x,y){x$group <- y; x}) %>%
  {function(x) do.call(rbind, x)}()

Y <- (DA_data$callow==1) + (DA_data$species=="F. sanguinea")

peak_prop <- peak_proportions_table(mass_spectra_data[mass_spectra_data$chromatogram_ID %in% DA_data$chromatogram_ID,])
peak_prop <- peak_prop[match(DA_data$chromatogram_ID, rownames(peak_prop)),]
# replace zeros with a small number before transformation
peak_prop[peak_prop==0] <- 10^-16
# apply central log ratio transformation
peak_transformed <- t(apply(peak_prop, 1, clr_transformation))
PCA_callow_discr <- prcomp(peak_transformed, scale. = TRUE)

load("R/sysdata.rda")
for(i in 1:20){
  if(!"true_sample_assignment_AUC" %in% ls()) {
    true_sample_assignment_AUC <- list()
    true_sample_assignment_AUC$AUC <- c()
    set.seed(3821)
  }
  else .Random.seed <- true_sample_assignment_AUC$random.seed
  discr_analysis_callow <- splsda(PCA_callow_discr$x, Y=Y ,multilevel=DA_data$group,
                                  ncomp=7, scale = FALSE)
  perf.discr_analysis_callow <- perf(discr_analysis_callow,
                                     validation = "Mfold",
                                     folds = 4,
                                     progressBar = FALSE,
                                     auc = TRUE,
                                     nrepeat = 10,
                                     scale = FALSE,
                                     cpus=6)
  discr_analysis_callows_tuned <- tune.splsda(PCA_callow_discr$x,
                                              Y=Y,
                                              multilevel=DA_data$group,
                                              ncomp = perf.discr_analysis_callow$choice.ncomp[1,1],
                                              test.keepX = seq(1:120),
                                              validation = 'Mfold',
                                              folds = 5, nrepeat = 15,
                                              dist = 'max.dist', # use max.dist measure
                                              measure = "BER",
                                              progressBar = FALSE,
                                              scale = FALSE,
                                              cpus=6)
  n_comp <- discr_analysis_callows_tuned$choice.ncomp$ncomp
  if(is.null(n_comp)) n_comp <- 1
  discr_analysis_callows_res <- splsda(PCA_callow_discr$x,Y=Y,
                                       multilevel= DA_data$group,
                                       keepX=discr_analysis_callows_tuned$choice.keepX,
                                       ncomp=n_comp,
                                       scale = FALSE)
  predicted_species <- calculate_SII(peak_prop, discr_analysis_callows_res, PCA_callow_discr, predicted_category=3)$predicted_species
  curve_data <- calculate_accuracy_metrics(predicted_species[Y>0], Y[Y>0]==2)
  AUC <- calculate_AUC(curve_data$FPR, curve_data$TPR)
  print(i)
  print(AUC)
  true_sample_assignment_AUC$AUC <- c(true_sample_assignment_AUC$AUC, AUC)
  true_sample_assignment_AUC$random.seed <- .Random.seed
  usethis::use_data(discriminant_analysis_randomization, true_sample_assignment_AUC, internal = TRUE, overwrite = TRUE)

}





for(i in 1:20){
  if(!"discriminant_analysis_randomization" %in% ls()) {
    discriminant_analysis_randomization <- list()
    discriminant_analysis_randomization$AUC <- c
    discriminant_analysis_randomization$Y <- Y
    set.seed(1924)
    print("Initialization")
  }
  else .Random.seed <- discriminant_analysis_randomization$random.seed
  Y <- discriminant_analysis_randomization$Y
  Y %>% split(.,DA_data$group) %>% lapply(function(x) {x[x>0]<-sample(x[x>0], length(x[x>0]));x}) %>% unlist() ->Y
  discr_analysis_callow <- splsda(PCA_callow_discr$x, Y=Y ,multilevel=DA_data$group,
                                  ncomp=7, scale = FALSE)
  perf.discr_analysis_callow <- perf(discr_analysis_callow,
                                     validation = "Mfold",
                                     folds = 4,
                                     progressBar = FALSE,
                                     auc = TRUE,
                                     nrepeat = 10,
                                     scale = FALSE,
                                     cpus=6)
  discr_analysis_callows_tuned <- tune.splsda(PCA_callow_discr$x,
                                              Y=Y,
                                              multilevel=DA_data$group,
                                              ncomp = perf.discr_analysis_callow$choice.ncomp[1,1],
                                              test.keepX = seq(1:120),
                                              validation = 'Mfold',
                                              folds = 5, nrepeat = 15,
                                              dist = 'max.dist', # use max.dist measure
                                              measure = "BER",
                                              progressBar = FALSE,
                                              scale = FALSE,
                                              cpus=6)
  n_comp <- discr_analysis_callows_tuned$choice.ncomp$ncomp
  if(is.null(n_comp)) n_comp <- 1
  discr_analysis_callows_res <- splsda(PCA_callow_discr$x,Y=Y,
                                       multilevel= DA_data$group,
                                       keepX=discr_analysis_callows_tuned$choice.keepX,
                                       ncomp=n_comp,
                                       scale = FALSE)
  predicted_species <- calculate_SII(peak_prop, discr_analysis_callows_res, PCA_callow_discr, predicted_category=3)$predicted_species
  curve_data <- calculate_accuracy_metrics(predicted_species[Y>0], Y[Y>0]==2)
  AUC <- calculate_AUC(curve_data$FPR, curve_data$TPR)
  print(i)
  print(AUC)
  discriminant_analysis_randomization$AUC <- c(discriminant_analysis_randomization$AUC, AUC)
  discriminant_analysis_randomization$random.seed <- .Random.seed
  discriminant_analysis_randomization$Y <- Y
  usethis::use_data(discriminant_analysis_randomization, true_sample_assignment_AUC, internal = TRUE, overwrite = TRUE)

}





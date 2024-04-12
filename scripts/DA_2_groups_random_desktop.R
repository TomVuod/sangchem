#setwd("/home/t.wlodarczyk/chemical_ecology/sangchem")
devtools::load_all(".")
args = 50
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

DA_data <- dplyr::filter(development_data, caste=="worker", !remarks %in% c("aggression actor", "aggression target", "queenless colony"),
                         species =="F. sanguinea") %>%
  dplyr::select(colony, species, callow, census_date, chromatogram_ID, sang_prop) %>%
  split(list(.$colony, .$census_date), drop=TRUE) %>%
  lapply(function(x){if(length(unique(x$callow))==1) return(NULL); x}) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  purrr::map2(seq_along(.), function(x,y){x$group <- y; x}) %>%
  {function(x) do.call(rbind, x)}()

Y <- as.numeric(DA_data$callow)

peak_prop <- peak_proportions_table(mass_spectra_data[mass_spectra_data$chromatogram_ID %in% DA_data$chromatogram_ID,])
peak_prop <- peak_prop[match(DA_data$chromatogram_ID, rownames(peak_prop)),]
# replace zeros with a small number before transformation
peak_prop[peak_prop==0] <- 10^-16
# apply central log ratio transformation
peak_transformed <- t(apply(peak_prop, 1, clr_transformation))
PCA_callow_discr <- prcomp(peak_transformed, scale. = TRUE)
PCA_callow_discr$x <- withinVariation(PCA_callow_discr$x, data.frame(DA_data$group))
results <- {if(file.exists("random_labels_2_groups.rds")) readRDS("random_labels_2_groups.rds")
  else {
    results <- list()
    results$AUC <- c()
    set.seed(1342)
    results$random.seed <- .Random.seed
    results$Y <- split(Y,DA_data$group) %>% lapply(function(x) {x<-sample(x, length(x));x}) %>% unlist() 
    results
  }
}
counter <- 0
while(TRUE){
  counter <- counter+1
  message(counter)
  .Random.seed <- results$random.seed
  Y <- results$Y
  discr_analysis_callow <- splsda(PCA_callow_discr$x, Y=Y,
                                  ncomp=7, scale = FALSE)
  perf.discr_analysis_callow <- perf(discr_analysis_callow,
                                     validation = "Mfold",
                                     folds = 4,
                                     progressBar = FALSE,
                                     auc = TRUE,
                                     nrepeat = 10,
                                     scale = FALSE,
                                     cpus=2)
  tryCatch(discr_analysis_callows_tuned <- tune.splsda(PCA_callow_discr$x,
                                                        Y=Y,
                                                        ncomp = perf.discr_analysis_callow$choice.ncomp[1,1],
                                                        test.keepX = seq(1:120),
                                                        validation = 'Mfold',
                                                        folds = 5, nrepeat = 15,
                                                        dist = 'max.dist', # use max.dist measure
                                                        measure = "BER",
                                                        progressBar = FALSE,
                                                        scale = FALSE,
                                                        cpus=2), error=function(cond) 0)

  n_comp <- discr_analysis_callows_tuned$choice.ncomp$ncomp
  if(is.null(n_comp)) n_comp <- 1
  discr_analysis_callows_res <- splsda(PCA_callow_discr$x,Y=Y,
                                       keepX=discr_analysis_callows_tuned$choice.keepX,
                                       ncomp=n_comp,
                                       scale = FALSE)
  predicted_species <- calculate_SII(peak_prop, discr_analysis_callows_res, PCA_callow_discr, predicted_category=1,
                                     multilevel=DA_data$group)$predicted_species
  curve_data <- calculate_accuracy_metrics(predicted_species, Y==0)
  AUC <- calculate_AUC(curve_data$FPR, curve_data$TPR)
  results$AUC <- c(results$AUC, AUC)
  print(AUC)
  results$random.seed <- .Random.seed
  results$Y <- split(Y,DA_data$group) %>% lapply(function(x) {x<-sample(x, length(x));x}) %>% unlist() 
  saveRDS(results, "random_labels_2_groups.rds")
  if(length(results$AUC)>=as.numeric(args[1])) break
}


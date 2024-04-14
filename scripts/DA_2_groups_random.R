setwd("/home/t.wlodarczyk/chemical_ecology/sangchem")
devtools::load_all(".")
args = commandArgs(trailingOnly=TRUE)
library(mixOmics)
library(dplyr)
# load data
data("development_data")
data("mass_spectra_data")

DA_data <- dplyr::filter(development_data, caste=="worker", !remarks %in% c("aggression actor", "aggression target", "queenless colony"),
                         species == "F. sanguinea") %>%
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
    set.seed(1342)
    results$random.seed <- .Random.seed
    results$Y <- list()
    results$predictions <- list()
    results
  }
}
counter=0
while(TRUE){
  counter <- counter+1
  message(counter)
if(length(results$Y)>0) Y <- results$Y[[length(results$Y)]]
  .Random.seed <- results$random.seed
  Y <- split(Y, DA_data$group) %>% lapply(function(x) {x<-sample(x, length(x));x}) %>% unlist()
  discr_analysis_callow <- splsda(PCA_callow_discr$x, Y=Y ,
                                  ncomp=7, scale = FALSE)
  perf.discr_analysis_callow <- perf(discr_analysis_callow,
                                     validation = "Mfold",
                                     folds = 4,
                                     progressBar = FALSE,
                                     auc = TRUE,
                                     nrepeat = 10,
                                     scale = FALSE,
                                     cpus=18)
  tryCatch({discr_analysis_callows_tuned <- tune.splsda(PCA_callow_discr$x,
                                                        Y=Y,
                                                        ncomp = perf.discr_analysis_callow$choice.ncomp[1,1],
                                                        test.keepX = seq(1:120),
                                                        validation = 'Mfold',
                                                        folds = 5, nrepeat = 15,
                                                        dist = 'max.dist', # use max.dist measure
                                                        measure = "BER",
                                                        progressBar = FALSE,
                                                        scale = FALSE,
                                                        cpus=18)

  n_comp <- discr_analysis_callows_tuned$choice.ncomp$ncomp
  if(is.null(n_comp)) n_comp <- 1
  discr_analysis_callows_res <- splsda(PCA_callow_discr$x,Y=Y,
                                       keepX=discr_analysis_callows_tuned$choice.keepX,
                                       ncomp=n_comp,
                                       scale = FALSE)

  perf.discr_analysis_callows_res <- perf(discr_analysis_callows_res,
                                     validation = "Mfold",
                                     folds = 4,
                                     progressBar = FALSE,
                                     nrepeat = 10,
                                     scale = FALSE,
                                     cpus=18)

  #curve_data <- calculate_accuracy_metrics(predicted_species[Y>0], Y[Y>0]==2)
  #AUC <- calculate_AUC(curve_data$FPR, curve_data$TPR)
  #results$AUC <- c(results$AUC, AUC)
  results$random.seed <- .Random.seed
  results$Y <- c(results$Y, list(Y))
  results$predictions <- c(results$predictions, perf.discr_analysis_callows_res$predict)
  saveRDS(results, "random_labels_2_groups.rds")}, error=function(cond) print(cond))
  if(length(results$predictions)>=as.numeric(args[1])) break
}


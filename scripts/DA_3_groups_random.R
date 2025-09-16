setwd("/home/t.wlodarczyk/chemical_ecology/sangchem")
source("./R/transformations.R")
source("./R/quality_control.R")
load("./data/development_data.rda")
load("./data/mass_spectra_data.rda")
args = commandArgs(trailingOnly=TRUE)
library(mixOmics)
library(dplyr)
# load data
data("development_data")
data("mass_spectra_data")
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

DA_data <- dplyr::filter(development_data, caste=="worker", !remarks %in% c("aggression actor", "aggression target", "queenless colony")) %>%
  dplyr::select(colony, species, callow, census_date, chromatogram_ID, sang_prop) %>%
  split(list(.$colony, .$census_date), drop=TRUE) %>%
  lapply(function(x){if(length(unique(x$callow))==1) return(NULL); x}) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  purrr::map2(seq_along(.), function(x,y){x$group <- y; x}) %>%
  {function(x) do.call(rbind, x)}()

Y <- as.numeric(DA_data$callow) + as.numeric(DA_data$species == "F. sanguinea")

peak_prop <- peak_proportions_table(mass_spectra_data[mass_spectra_data$chromatogram_ID %in% DA_data$chromatogram_ID,])
peak_prop <- peak_prop[match(DA_data$chromatogram_ID, rownames(peak_prop)),]
peak_prop <- peak_prop[,apply(peak_prop,2,function(x) sum(x)>0)]
peak_prop <- freq_abun_QC(peak_prop, freq_cutoff = 0.3, abundance_cutoff = 0.005)
# replace zeros with a small number before transformation
peak_prop[peak_prop==0] <- 10^-5
# apply central log ratio transformation
peak_transformed <- t(apply(peak_prop, 1, clr_transformation))
PCA_callow_discr <- prcomp(peak_transformed, scale. = TRUE)
n_PCs <- sum(cumsum(apply(PCA_callow_discr$x,2,var))/sum(apply(PCA_callow_discr$x,2,var))<0.8)+1
PCA_callow_discr$rotation <- PCA_callow_discr$rotation[,1:n_PCs]
PCA_callow_discr$x <- PCA_callow_discr$x[,1:n_PCs]
file_name <- paste0("./output/random_labels_3_groups_",as.character(args[1]),".rds")
results <- {if(file.exists(file_name)) readRDS(file_name)
  else {
    results <- list()
    set.seed(1342)
    rseed <- round(runif(10)*1000)[as.numeric(args[1])]
    set.seed(rseed)
    results$random.seed <- .Random.seed
    results$Y <- Y
    results$predictions <- list()
    results
  }
}
counter=0
while(TRUE){
  counter <- counter+1
  message(counter)
  Y <- results$Y
  .Random.seed <- results$random.seed
  Y <- split(Y, DA_data$group) %>% lapply(function(x) {x[x>0]<-sample(x[x>0], length(x[x>0]));x}) %>% unlist()
  discr_analysis_callow <- splsda(PCA_callow_discr$x, Y=Y, multilevel = DA_data$group,
                                  ncomp=7, scale = FALSE)
  tryCatch({perf.discr_analysis_callow <- perf(discr_analysis_callow,
                                               validation = "Mfold",
                                               folds = 4,
                                               progressBar = FALSE,
                                               auc = TRUE,
                                               nrepeat = 30,
                                               scale = FALSE,
                                               BPPARAM = BiocParallel::MulticoreParam(RNG=sample(9999)))
  selected_n <- perf.discr_analysis_callow$choice.ncom["overall", "max.dist"]
  discr_analysis_callows_tuned <- tune.splsda(PCA_callow_discr$x,
                                              Y=Y,
                                              ncomp = selected_n,
                                              test.keepX = seq(1:dim(PCA_callow_discr$x)[2]),
                                              validation = 'Mfold',
                                              folds = 4, nrepeat = 30,
                                              multilevel = DA_data$group,
                                              dist = 'max.dist', # use max.dist measure
                                              measure = "BER",
                                              progressBar = FALSE,
                                              scale = FALSE,
                                              BPPARAM = BiocParallel::MulticoreParam(RNG=sample(9999)))

  selected_n <- discr_analysis_callows_tuned$choice.ncomp$ncomp
  if(is.null(n_comp)) n_comp <- 1
  discr_analysis_callows_res <- splsda(PCA_callow_discr$x,Y=Y,
                                       keepX=discr_analysis_callows_tuned$choice.keepX[1:selected_n],
                                       ncomp=selected_n,multilevel = DA_data$group,
                                       scale = FALSE)

  perf.discr_analysis_callows_res <- perf(discr_analysis_callows_res,
                                          validation = "Mfold",
                                          folds = 4,
                                          progressBar = FALSE,
                                          nrepeat = 1,
                                          scale = FALSE,
                                          BPPARAM = BiocParallel::MulticoreParam(RNG=sample(9999)))

  curve_data <- calculate_accuracy_metrics(perf.discr_analysis_callows_res$predict[[length(perf.discr_analysis_callows_res$predict)]][,3,1][Y>0], Y[Y>0]==2)
  AUC <- calculate_AUC(curve_data$FPR, curve_data$TPR)
  results$AUC <- c(results$AUC, AUC)
  results$random.seed <- .Random.seed
  results$Y <- Y
  saveRDS(results, file_name)}, error=function(cond) message(cond))
  if(length(results$AUC)>=as.numeric(100)) break
}


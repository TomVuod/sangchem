##############
# distances_to_free_living_fusca
##############

set.seed(10101)
# Calculate emprical p-value for the treatment 1
res_1 <- distances_permutation_test(treatment_id=1, exclude_naked_pupae=TRUE)

# Calculate emprical p-value for the treatment 2
res_2 <- distances_permutation_test(treatment_id=2, exclude_naked_pupae=TRUE)

# Calculate emprical p-value for the treatment 3
res_3 <- distances_permutation_test(treatment_id=3, exclude_naked_pupae=TRUE)

# Calculate emprical p-value for the treatment 3
res_4 <- distances_permutation_test(treatment_id=4, exclude_naked_pupae=TRUE)

distances_to_free_living_fusca <- list()
distances_to_free_living_fusca$cocoons <- list()
distances_to_free_living_fusca$cocoons$treatment_1 <- res_1
distances_to_free_living_fusca$cocoons$treatment_2 <- res_2
distances_to_free_living_fusca$cocoons$treatment_3 <- res_3
distances_to_free_living_fusca$cocoons$treatment_4 <- res_4

set.seed(10101)
# Calculate emprical p-value for the treatment 1
res_1 <- distances_permutation_test(treatment_id=1)

# Calculate emprical p-value for the treatment 2
res_2 <- distances_permutation_test(treatment_id=2)

# Calculate emprical p-value for the treatment 3
res_3 <- distances_permutation_test(treatment_id=3)

# Calculate emprical p-value for the treatment 3
res_4 <- distances_permutation_test(treatment_id=4)

distances_to_free_living_fusca$all<- list()
distances_to_free_living_fusca$all$treatment_1 <- res_1
distances_to_free_living_fusca$all$treatment_2 <- res_2
distances_to_free_living_fusca$all$treatment_3 <- res_3
distances_to_free_living_fusca$all$treatment_4 <- res_4
usethis::use_data(distances_to_free_living_fusca, overwrite = TRUE)

##############
# globals
##############

suppressMessages({
  library(stringi)
  library(dplyr)
  library(mixOmics)
})

devtools::load_all()

# load data from the the earlier study
data("field_colonies")

# select samples from F. sanguinea ants and free-living F. fusca ants encoded in chromatogram_ID
fusca_samples <- as.character(unique(field_colonies$chromatogram_ID)[stri_detect_regex(unique(field_colonies$chromatogram_ID), "fus[0-9].?-[0-9]\\([1-2]\\)")])
sanguinea_samples<-as.character(unique(field_colonies$chromatogram_ID)[stri_detect_regex(unique(field_colonies$chromatogram_ID),"sangFs[0-9].?\\([0-9]\\)")])

# remove outliers
sanguinea_samples <- setdiff(sanguinea_samples,c("sangFs2(1)","sangFs5(2)","sangFs26(1)","sangFs26(2)"))
fusca_samples <- setdiff(fusca_samples,c("fus9-2(1)", "fus26-3(2)"))
discrimination_data <- filter(field_colonies,chromatogram_ID %in% c(fusca_samples,sanguinea_samples))
peak_prop_unfiltered <- peak_proportions_table(discrimination_data)
# feature selection
peak_prop <- freq_abun_QC(peak_prop_unfiltered)
# re-normalize rows to 1
peak_prop <- peak_prop/rowSums(peak_prop)
# replace zeros with a small number before transformation
peak_prop[peak_prop==0] <- 10^-16
# apply central log ratio transformation
peak_transformed <- t(apply(peak_prop, 1, clr_transformation))

PCA_species_discr <- prcomp(peak_transformed, scale. = TRUE)

Y <- rownames(PCA_species_discr$x) %in% sanguinea_samples

set.seed(9230)
# run the main model
species_discrimination_PC.plsda <- plsda(PCA_species_discr$x, Y, ncomp=10, scale=FALSE)
plotIndiv(species_discrimination_PC.plsda,group=Y,
          legend=TRUE,ellipse=TRUE,comp=1:2)

# determine the optimal number of components used in the model
spec_discr_PC.plsda.perf<-perf(species_discrimination_PC.plsda, nrepeat = 10,folds=5,
                               progressBar = FALSE,validation = "Mfold",
                               auc=TRUE, scale=FALSE)


selected_n <- spec_discr_PC.plsda.perf$choice.ncom["overall", "max.dist"]

# determine the number of the original variables to be kept in each component (actually
# here we have only one component)
spec_discr_PC.splsda <- tune.splsda(PCA_species_discr$x,Y,test.keepX = c(1:dim(PCA_species_discr$x)[2]),ncomp = selected_n,
                                    validation="Mfold", progressBar=FALSE, dist="max.dist",
                                    nrepeat=5, cpus=2, folds=5, measure="BER", scale=FALSE)

# final model
species_prediction_model <- splsda(PCA_species_discr$x, Y, ncomp=selected_n,
                                   keepX =spec_discr_PC.splsda$choice.keepX, scale=FALSE)

# visualize results for the first two components
plotIndiv(species_prediction_model, group=Y, legend=TRUE, ellipse=TRUE, comp=1:2)

# quantify predicted species (F. sanguinea closer to one and F. fusca closer to 0)
predicted_species <- predict(species_prediction_model, PCA_species_discr$x)$predict[,2,species_prediction_model$ncomp]

library(Hmisc)
cor_p_val<-rcorr(cbind(peak_prop, predicted_species))
cor_p_val$P[cor_p_val$r==1]<-0
#significance threshold = 0.05
significant_correlations <- cor_p_val$r["predicted_species",][cor_p_val$P["predicted_species",]<0.05]
peaks_sang <- names(significant_correlations)[significant_correlations>0]
peaks_sang <- as.numeric((setdiff(peaks_sang,"predicted_species")))
peaks_fusca <- as.numeric(names(significant_correlations)[significant_correlations<0])
print(peaks_sang)
print(peaks_fusca)

globals <- readRDS(system.file("globals.rds", package="sangchem"))
globals[["species_prediction_model"]] <- species_prediction_model
globals[["PCA_species_discr"]] <- PCA_species_discr
globals[["peak_sang"]] <- peaks_sang
globals[["peak_fusca"]] <- peaks_fusca
saveRDS(globals, system.file("globals.rds", package="sangchem"))

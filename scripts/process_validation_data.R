# process cross-validation on high-performance computing nodes
# merge files with partial results into single ones

setwd("/home/t.wlodarczyk/chemical_ecology/sangchem")
AUC_basic <- c()
for(i in 1:10) AUC_basic <- c(AUC_basic, readRDS(paste0("./output/true_labels_3_groups_", i, ".rds"))$AUC[1:100])

AUC_random <- c()
for(i in 1:10) AUC_random <- c(AUC_random, readRDS(paste0("./output/random_labels_3_groups_", i, ".rds"))$AUC[1:100])


discrimination_model_validation <- list(AUC_observed = AUC_basic, AUC_null = AUC_random)
saveRDS(discrimination_model_validation, "./data/discrimination_model_validation.rda")

CHC_mass_prediction <- readRDS("./output/predicted_values.rds")
saveRDS(CHC_mass_prediction, "./data/CHC_mass_prediction.rda")

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


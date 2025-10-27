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

annotate_peaks<-function(data, ref, parms=data.frame(tolerance=c(0.06, 0.06, 0.045, 0.03,rep(c(0.05, 0.04, 0.03), 100))[1:100], annotation_scope=rep(1.5, 100)[1:100],
                                                     calculation_scope=rep(2, 100)[1:100], step=c(0,0,0,1,rep(c(0,0,1),100))[1:100]),
                         start_time=32, end_time=58, accuracy=0.02, fixed_peaks = c(), ...){
  annotation_scope <- 0
  pointA <- start_time
  counter <- 0
  temp_data <- data.frame()
  while(pointA + annotation_scope < end_time){
    if (counter > 0){
      temp_data <- temp_data[(temp_data[,"retention_time"]>=(pointA-annotation_scope))&(temp_data[,"retention_time"]<=(pointA+annotation_scope)),]
      data <- rbind(data,temp_data)
      pointA <- pointA + step
    }
    print(sprintf("Punkt środkowy: %1.2f", pointA))
    counter<-counter+1
    if (counter>nrow(parms)) stop("parms wyczerpane")
    tolerance <- parms$tolerance[counter]
    step <- parms$step[counter]
    annotation_scope <- parms$annotation_scope[counter]
    calculation_scope <- parms$calculation_scope[counter]
    if (annotation_scope>calculation_scope) warning('Zasieg obliczen mniejszy od zasiegu oznaczania')
    temp_data <- data[(data[,"retention_time"]>=(pointA-calculation_scope))&(data[,"retention_time"]<=(pointA+calculation_scope)),]
    data <- data[(data[,"retention_time"]<(pointA-annotation_scope))|(data[,"retention_time"]>(pointA+annotation_scope)),]
    shared_peaks <- base::intersect(ref[,"peak_ID"], temp_data[,"peak_ID"])
    shared_peaks <- shared_peaks[!is.na(shared_peaks)&shared_peaks != 0]
    print(do.call(sprintf,c(list(paste0("Liczba wspólnych pików: %d. Piki: ", paste(rep(" %d, ", length(shared_peaks)),collapse='')),length(shared_peaks)),as.list(shared_peaks))))
    if (length(shared_peaks)<2&length(shared_peaks)>0) {
      czas1<-calculate_mean_rt(ref, peak_subset=shared_peaks)
      czas2<-calculate_mean_rt(temp_data, peak_subset=shared_peaks)
      czasy_tab<-calculate_mean_rt(ref, add_peaks = TRUE)
      czasy_tab<-cbind(czasy_tab, czasy_przewidywane=((czasy_tab[,"retention_time"])*rep((czas2/czas1), nrow(czasy_tab))))
      temp_data<-temp_data[(temp_data[,"retention_time"]>=(pointA-annotation_scope))&(temp_data[,"retention_time"]<=(pointA+annotation_scope)),]
      if (nrow(temp_data)==0) next
      temp_data <- recognize_peaks(temp_data, czasy_tab, tolerance=tolerance,
                                   fixed_peaks_ = fixed_peaks)
      next
    }
    if (length(shared_peaks)==0) next
    print("Wywołanie regresji czasów")
    res<-predict_ret_times(temp_data, ref,
                           max_deviation=accuracy, fixed_peaks = fixed_peaks, ...)
    print("Czasy przewidywane")
    print(res)
    print("Czasy oznaczonych pików")
    print(calculate_mean_rt(temp_data, add_peaks = TRUE))
    if (is.null(res)) next
    temp_data <- temp_data[(temp_data[,"retention_time"]>=(pointA-annotation_scope))&(temp_data[,"retention_time"]<=(pointA+annotation_scope)),]
    peaks_split <- setdiff(intersect(temp_data$peak_ID, data$peak_ID), 0)
    data <- rbind(data, temp_data[temp_data$peak_ID %in% peaks_split,])
    temp_data <- temp_data[!temp_data$peak_ID %in% peaks_split,]
    if (nrow(temp_data)==0) next
    temp_data <- recognize_peaks(temp_data, res, tolerance=tolerance,
                                 fixed_peaks_ = fixed_peaks,
                                 excluded_peaks = peaks_split)
    print("Piki po oznaczeniu:")
    print(unique(temp_data$peak_ID))
    if (sum(temp_data[,"peak_ID"]>0)>0) { #temp_data są wprowadzana jako argument "referencje" do funkcji
      #"numeracja pików" dlatego powinny zawierać przynajmniej jeden oznaczony pik
      temp_data <- recognize_peaks(temp_data, temp_data, tolerance=tolerance, fixed_peaks_ = fixed_peaks)
    }
  }
  temp_data<-temp_data[(temp_data[,"retention_time"]>=(pointA-annotation_scope))&(temp_data[,"retention_time"]<=(pointA+annotation_scope)),]
  rbind(data,temp_data)
}

lm_2 <- function(x,y,weights,squared_term = FALSE){
  if (squared_term)
    return(lm(y~x + I(x^2), weights = weights))
  return(lm(y~x, weights = weights))
}

pred_2 <- function(model, x){
  coeffs <- coef(model)
  if(length(coeffs) == 3)
    return(coeffs[1]+coeffs[2]*x + coeffs[3]*x^2)
  else
    coeffs[1]+coeffs[2]*x
}

recognize_peaks <- function(data, referencje, tolerance=0.03, fixed_peaks_ = c(),
                            excluded_peaks = NULL,
                            z.max = 1, ...) {
  if (!is.data.frame(data)) stop("Argument is not a data frame obejct")
  if (nrow(data) == 0) stop("Empty data frame")
  protected_data <- data[(data$peak_ID %in% fixed_peaks_), ]
  data <- data[!(data$peak_ID %in% fixed_peaks_), ]
  if (!is.data.frame(referencje)) stop("Argument is not a data frame obejct")
  if (nrow(data)==0) return(rbind(data, protected_data))
  if (nrow(data)==0 & nrow(protected_data)==0) stop("Numeracja pików: brak danych")
  referencje <- referencje[referencje$peak_ID!=0,]
  mean_rts <- calculate_mean_rt(referencje, add_peaks=TRUE)
  if (nrow(mean_rts)==0) return(rbind(data, protected_data))
  mean_rts <- mean_rts[!mean_rts$peak_ID %in% fixed_peaks_,]
  if(nrow(mean_rts)==0) return(rbind(data, protected_data))
  data$peak_ID <- 0
  for(i in 1:nrow(data)){
    time_diff = abs(data[i,"retention_time"]-mean_rts$retention_time)
    min_diff_ind = which.min(time_diff)
    if (time_diff[min_diff_ind]>tolerance) next
    if (data[i,"peak_ID"] %in% c(fixed_peaks_, excluded_peaks)) next
    data[i, "peak_ID"] <- mean_rts$peak_ID[min_diff_ind]
  }
  mutable_peaks <- setdiff(data$peak_ID, c(fixed_peaks_, 0))
  n_sd <- list(...)$n_sd
  if(is.null(n_sd)) n_sd = 2
  data <- refine_peaks(data, mutable_peaks, n_sd = n_sd)
  # for(peak in mutable_peaks){
  #   rts <- data[data$peak_ID==peak, "retention_time"]
  #   if(length(rts)<2) next
  #   z <- abs(rts-mean(rts))/sd(rts)
  #   data[data$peak_ID==peak,]$peak_ID[z>z.max] <- 0
  #   print("z max")
  #   rts <- data[data$peak_ID==peak, "retention_time"]
  #   print(max(abs(rts-mean(rts))/sd(rts)))
    # while(TRUE){
    #   rts <- data[data$peak_ID==peak, "retention_time"]
    #   if(length(rts)<2) break
    #   z <- abs(rts-mean(rts))/sd(rts)
    #   if(max(z)<z.max) break
    #   data[data$peak_ID==peak,]$peak_ID[which.max(x)] <- 0
    # }
  # }
  rbind(data, protected_data)
}

predict_ret_times <- function(data, reference, max_deviation = 0.04,
                              fixed_peaks = c(), squared_term = FALSE){
  library(dplyr)
  reference <- filter(reference, peak_ID > 0)
  if (any(is.na(c(data$retention_time, reference$retention_time)))) stop("NA values for retention time")
  shared_peaks<-sort(base::intersect(data$peak_ID, reference$peak_ID))
  shared_peaks<-shared_peaks[shared_peaks!=0&!is.na(shared_peaks)]
  data_times<-calculate_mean_rt(data, peak_subset = shared_peaks)
  ref_times<-calculate_mean_rt(reference, peak_subset=shared_peaks, add_peaks = TRUE)
  peaks_n <- table(reference$peak_ID)
  if(length(shared_peaks) == 2) squared_term <- FALSE
  while(TRUE){
    lmodel <- lm_2(ref_times$retention_time, data_times, weights = peaks_n[as.character(ref_times$peak_ID)], squared_term)
    print(data.frame(peak = shared_peaks, predicted = lmodel$fitted.values,
                     resid = abs(lmodel$residuals)))
    if (all(shared_peaks %in% fixed_peaks)) break
    if (sum(abs(lmodel$residuals)[!(shared_peaks %in% fixed_peaks)] > max_deviation) == 0) break
    resid_decreasing = order(abs(lmodel$residuals), decreasing = TRUE)
    max_val_ind <- resid_decreasing[!shared_peaks[resid_decreasing] %in% fixed_peaks][1]
    data_times <- data_times[-max_val_ind]
    ref_times <- ref_times[-max_val_ind,]
    shared_peaks<-shared_peaks[-max_val_ind]
    if (length(shared_peaks)<(2+squared_term)) {
      print("To few peaks for regression")
      return(NULL)
    }
  }
  mean_rt <- calculate_mean_rt(reference, add_peaks = TRUE)
  pred_times <- pred_2(lmodel, mean_rt$retention_time)
  res <- data.frame(retention_time = pred_times, peak_ID = mean_rt$peak_ID)
  res[order(res$retention_time),]
}


calculate_mean_rt<-function(data, add_peaks = FALSE, peak_subset = NULL) {
  library(dplyr)
  if (!is.data.frame(data)) stop("Object is not a data frame")
  if(nrow(data)==0) stop("Empty data frame")
  # apply peak_ID order as given in peak_subset
  if (!is.null(peak_subset)) data <- filter(data, peak_ID %in% peak_subset)
  res <- group_by(data, peak_ID) %>%
    summarise(retention_time = mean(retention_time)) %>% arrange(peak_ID) %>%
    as.data.frame()
  if (!add_peaks) res <- pull(res, retention_time)
  res
}

library(spatstat)
refine_peaks <- function(data_df, selected_peaks = NULL, n_sd = 2){
  if(is.null(selected_peaks)) selected_peaks <- setdiff(data_df$peak_ID, 0)
  if(length(selected_peaks)==0) return(data_df)
  for(peak in selected_peaks){
    rts <- data_df[data_df$peak_ID == peak, "retention_time"]
    if(length(rts)<3) next
    rt_points <- ppp(rts, rep(0, length(rts)), range(rts), c(0,0))
    distances <- pairdist(rt_points)
    dist_threshold <- sd(rts)*n_sd
    adj_matrix <- pairdist(rt_points) <= dist_threshold
    diag(adj_matrix) <- FALSE
    gr <- igraph::graph_from_adjacency_matrix(adj_matrix)
    ind <- as.numeric(igraph:::largest_cliques(gr)[[1]])
    data_df[data_df$peak_ID == peak,]$peak_ID[-ind] <- 0
  }
  data_df
}

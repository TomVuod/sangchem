# peak order
batch7 <- mass_spectra_data[mass_spectra_data$batch==7 & mass_spectra_data$peak_ID!=0,]
setdiff(mass_spectra_data$peak_ID , batch7$peak_ID)
peaks <- unique(batch7$peak_ID)
mean_ret_times <- unlist(lapply(peaks, function(x) mean(batch7[batch7$peak_ID==x, "retention_time"])))
peaks_ordered <- peaks[order(mean_ret_times)]

peak_ID_df <- data.frame(peak_ID = 1:length(peaks_ordered), old_ID = peaks_ordered)
peak_ID_df$identification <- unlist(lapply(peaks_ordered, function(x) mass_spectra_data[mass_spectra_data$peak_ID==x,"identification"][1]))

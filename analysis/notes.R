# change to development_data

development_data[development_data$colony=="SD18-11"&development_data$sang_prop>0.1,"remarks"] <- "queenless colony"


suppressMessages({
  library(sangchem)
  library(stringi)
  library(dplyr)
  library(mixOmics)
})

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

unit_input <- diag(rep(1, ncol(peak_prop)))
unit_input[unit_input==0] <- 10^-16
colnames(unit_input) <- rownames(unit_input) <-colnames(peak_prop)

res <- calculate_SII(unit_input)



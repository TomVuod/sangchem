# choose first 10 samples per species to analyse, selecting them so that to 
# evently subsample the empirical distribution of sizes
library(sangchem)
data("development_data")
head_widths <- development_data[,c("chromatogram_ID", "species", "head_width")]
# F. sanguinea
head_widths_sang <- head_widths[head_widths$species == "F. sanguinea",]
head_widths_sang <- head_widths_sang[!is.na(head_widths_sang$head_width),]
# determine indices of ordered samples to be chosen
indices <- findInterval(1:10, seq(1, 10, length.out = nrow(head_widths_sang)))
chr_IDs_sang <- head_widths_sang[order(head_widths_sang$head_width),][indices,"chromatogram_ID"]

# F. fusca
head_widths_fus <- head_widths[head_widths$species == "F. fusca",]
head_widths_fus <- head_widths_fus[!is.na(head_widths_fus$head_width),]
# determine indices of ordered samples to be chosen
indices <- findInterval(1:10, seq(1, 10, length.out = nrow(head_widths_fus)))
chr_IDs_fus <- head_widths_fus[order(head_widths_fus$head_width),][indices,"chromatogram_ID"]

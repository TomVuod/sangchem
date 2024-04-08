hw_data <- read.csv("C:/Users/uwb/Documents/Badania/Ekologia chemiczna kolonii inicjalnych F. sanguinea/Data analysis/Raw data/head_width.csv")


all_samples <- sort(unique(hw_data$chromatogram_ID))


for(i in seq_len(floor(length(all_samples)/16))){
  index_first = (i-1)*16+1
  index_last = min((i*16), length(all_samples))
  cat(sprintf("Kubek SzG %d: Eppendrof od %d do %d\n", i, all_samples[index_first], all_samples[index_last]))
}

cat("Reszta do kubka SzG 0")
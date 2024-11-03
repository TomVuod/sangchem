devtools::load_all()
data("development_data")
data("separation_experiment_data")
res_table <- read.csv(system.file("Data-raw/surface_area.csv",package="sangchem"))
head_widths <- c()
species <- c()
callow <- c()
for(i in 1:nrow(res_table)){
  chr_ID <- res_table$ID[i]
  if(chr_ID %in% development_data$chromatogram_ID){
    head_widths[i] <- as.numeric(development_data[development_data$chromatogram_ID==chr_ID,"head_width"])
    species[i] <- development_data[development_data$chromatogram_ID==chr_ID,"species"]
    callow[i] <- development_data[development_data$chromatogram_ID==chr_ID,"callow"]
    }
  else {
    head_widths[i] <- as.numeric(separation_experiment_data[separation_experiment_data$chromatogram_ID==chr_ID,"head_width"])
    species[i] <- "F. sanguinea"
    callow[i] <- NA
  }
}

res_table$species <- species
res_table$head_width <- head_widths
res_table[res_table$ID==389,"head_width"] = 946 # measured again
res_table$normalized_area <- res_table$area/head_widths^2
res_table$callow <- callow
res_table <- res_table[res_table$ID!=356,]
# sample 356 should be F. fusca but is F. sanguinea; maybe it is sample 350
# measurement of sample 31_h should be repeated; the same about 254_t


library(ggplot2)
ggplot(res_table, aes(y=normalized_area, x=view, fill=species))+
  geom_boxplot()
res_table |> split(list(res_table$species, res_table$view)) |> lapply(nrow)
mean_areas <- res_table |> split(list(res_table$species, res_table$view)) |> lapply(function(x) mean(x$normalized_area))
mean_areas[[1]]/mean_areas[[2]]
mean_areas[[3]]/mean_areas[[4]]
mean_areas[[5]]/mean_areas[[6]]

library(BiocParallel)
library(dplyr)
library(lme4)
setwd("/home/t.wlodarczyk/chemical_ecology/sangchem")
load("./data/development_data.rda")
load("./data/mass_spectra_data.rda")
load("./data/body_surface_area.rda")
source("./R/body_area_normalization.R")
CHC_normalizer <- get_body_area_normalizer(body_surface_area, development_data)
dev_samples <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor", "aggression target", "queenless colony"), callow==0) %>%
  dplyr::select(chromatogram_ID, species, colony, sang_prop, head_width, mass, census_date) %>%
  mutate(normalized_mass=mass*CHC_normalizer(species, head_width))

model_input <- dev_samples[dev_samples$species=="F. sanguinea",]
# remove outliers
model_input <- model_input[-which.max(model_input$normalized_mass),]
lm_model <- lmer(I(log(normalized_mass)) ~ sang_prop   + (1|colony) +(1 | colony:census_date), data=model_input)

simulation_number=10^3

predict_values <- function(n, model, sang_prop){
  exp(bootMer(model,function(x) predict(x,data.frame(sang_prop=sang_prop),re.form=NA),nsim=simulation_number)$t)
}
pred_values_sang <- unlist(bplapply(1:1000, predict_values, model = lm_model, sang_prop = 1, BPPARAM = MulticoreParam(RNG=4381)))

model_input <- dev_samples[dev_samples$species=="F. fusca",]
lm_model <- lmer(I(log(normalized_mass)) ~ sang_prop  + (1|colony), data=model_input)
pred_values_fusca <- unlist(bplapply(1:1000, predict_values, model = lm_model, sang_prop = 0, BPPARAM = MulticoreParam(RNG=4381)))

res = list(pred_values_sang=pred_values_sang, pred_values_fusca=pred_values_fusca)
save(res, file= "/home/t.wlodarczyk/chemical_ecology/sangchem/data/CHC_mass_prediction.rda")


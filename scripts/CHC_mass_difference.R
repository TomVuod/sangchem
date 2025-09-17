library(BiocParallel)
library(dplyr)
library(lme4)
setwd("/home/t.wlodarczyk/chemical_ecology/sangchem")
load("./data/development_data.rda")
load("./data/mass_spectra_data.rda")
dev_samples <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor", "aggression target", "queenless colony"), callow==0) %>%
  dplyr::select(chromatogram_ID, species, colony, sang_prop, head_width, mass, census_date) %>%
  mutate(corrected_mass=mass/head_width^2*10^6)

model_input <- dev_samples[dev_samples$species=="F. sanguinea",]
# remove outliers
model_input <- model_input[-which.max(model_input$corrected_mass),]
lm_model <- lmer(I(log(corrected_mass)) ~ sang_prop  + I(head_width^2) + (1|colony) +(1 | colony:census_date), data=model_input)

simulation_number=10^3

predict_values <- function(n, model, sang_prop){
  exp(bootMer(model,function(x) predict(x,data.frame(sang_prop=sang_prop, head_width=1300),re.form=NA),nsim=simulation_number)$t)
}
pred_values_sang <- unlist(bplapply(1:1000, predict_values, model = lm_model, sang_prop = 1, BPPARAM = MulticoreParam(RNG=4381)))

model_input <- dev_samples[dev_samples$species=="F. fusca",]
lm_model <- lmer(I(log(corrected_mass)) ~ sang_prop + I(head_width^2)  + (1|colony), data=model_input)
pred_values_fusca <- unlist(bplapply(1:1000, predict_values, model = lm_model, sang_prop = 0, BPPARAM = MulticoreParam(RNG=4381)))


saveRDS(list(pred_values_sang=pred_values_sang, pred_values_fusca=pred_values_fusca), "/home/t.wlodarczyk/chemical_ecology/sangchem/output/predicted_values.rds")


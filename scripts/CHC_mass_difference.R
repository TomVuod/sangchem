  library(sangchem)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(DHARMa)
  library(BiocParallel)
data("mass_spectra_data")
data("development_data")
dev_samples <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor", "aggression target", "queenless colony"), callow==0) %>%
  dplyr::select(chromatogram_ID, species, colony, sang_prop, head_width, mass, census_date) %>%
  mutate(corrected_mass=mass/head_width^2*10^6)

model_input <- dev_samples[dev_samples$species=="F. sanguinea",]
lm_model <- lmer(I(log(corrected_mass)) ~ sang_prop  + (1|colony) +(1 | colony:census_date), data=model_input)


simulation_number=10^6
seeds <- (1:1000)*17
predict_values <- function(n, model, sang_prop, seeds){
  set.seed(seeds[n])
  exp(bootMer(model,function(x) predict(x,data.frame(sang_prop=sang_prop),re.form=NA),nsim=10^3)$t)
}
pred_values_sang <- bplapply(1:10000, predict_values, model = lm_model, sang_prop = 1, seeds = seeds, BPPARAM = MulticoreParam())
pred_values_sang <- do.call(c, pred_values_sang)

model_input <- dev_samples[dev_samples$species=="F. fusca",]
lm_model <- lmer(I(log(corrected_mass)) ~ sang_prop  + (1|colony), data=model_input)
pred_values_fusca <- bplapply(1:10000, predict_values, model = lm_model, sang_prop = 0, seeds = seeds, BPPARAM = MulticoreParam())
pred_values_fusca <- do.call(c, pred_values_fusca)

saveRDS(list(pred_values_sang=pred_values_sang, pred_values_fusca=pred_values_fusca), "/home/t.wlodarczyk/chemical_ecology/sangchem/predicted_values.rds")


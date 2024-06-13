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
  rseeds <- list()
  for(seed_int in (1:1000)*17){
    set.seed(seed_int)
    rseeds <- c(rseeds, list(.Random.seed))
  }
  seeds_fusca <- seeds_sang <- rseeds
if(file.exists("/home/t.wlodarczyk/chemical_ecology/sangchem/predicted_values.rds")){
  predicted_values <- readRDS("/home/t.wlodarczyk/chemical_ecology/sangchem/predicted_values.rds")
  pred_values_sang <- predicted_values$pred_values_sang
  pred_values_fusca <- predicted_values$pred_values_fusca
  seeds_sang <- seeds_sang
  seeds_fusca <- seeds_fusca
}
model_input <- dev_samples[dev_samples$species=="F. sanguinea",]
lm_model <- lmer(I(log(corrected_mass)) ~ sang_prop  + (1|colony) +(1 | colony:census_date), data=model_input)


simulation_number=10^4

seeds <- (1:1000)*17
predict_values <- function(n, model, sang_prop, seeds){
  set.seed(seeds[n])
  list(t=exp(bootMer(model,function(x) predict(x,data.frame(sang_prop=sang_prop),re.form=NA),nsim=simulation_number)$t),
       rseed = .Random.seed)
}
pred_values_sang <- bplapply(1:1000, predict_values, model = lm_model, sang_prop = 1, seeds = seeds, BPPARAM = MulticoreParam())
res_sang <- c()
seeds_sang <- list()
for(i in seq_along(pred_values_sang)){
  res_sang <- c(res_sang, pred_values_sang[[i]]$t)
  seeds_sang <- c(seeds_sang, pred_values_sang[[i]][2])
}

model_input <- dev_samples[dev_samples$species=="F. fusca",]
lm_model <- lmer(I(log(corrected_mass)) ~ sang_prop  + (1|colony), data=model_input)
pred_values_fusca <- bplapply(1:1000, predict_values, model = lm_model, sang_prop = 0, seeds = seeds, BPPARAM = MulticoreParam())
res_fusca <- c()
seeds_fusca <- list()
for(i in seq_along(pred_values_fusca)){
  res_fusca <- c(res_fusca, pred_values_fusca[[i]]$t)
  seeds_fusca <- c(seeds_fusca, pred_values_fusca[[i]][2])
}

saveRDS(list(pred_values_sang=res_sang, pred_values_fusca=res_fusca, seeds_sang = seeds_sang, seeds_fusca = seeds_fusca), "/home/t.wlodarczyk/chemical_ecology/sangchem/predicted_values.rds")


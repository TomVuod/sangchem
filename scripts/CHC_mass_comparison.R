suppressMessages({
  library(sangchem)
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(DHARMa)
})
# load data
data("development_data")
load_globals()
dev_samples <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor", "aggression target", "queenless colony"), callow==0) %>%
  dplyr::select(chromatogram_ID, species, colony, sang_prop, head_width, mass, census_date) %>%
  mutate(corrected_mass=mass/head_width^2*10^6)

model_input <- dev_samples[dev_samples$species=="F. sanguinea",]
lm_model <- lmer(I(log(corrected_mass)) ~ sang_prop  + (1|colony) +(1 | colony:census_date), data=model_input)
print(summary(lm_model))

simulation_number=10^3

sang_fam_simulated<-bootMer(lm_model,function(x) predict(x,data.frame(sang_prop=1),re.form=NA),nsim=simulation_number)$t
pred_values_sang <- exp(sang_fam_simulated)


model_input <- dev_samples[dev_samples$species=="F. fusca",]
lm_model <- lmer(I(log(corrected_mass)) ~ sang_prop  + (1|colony), data=model_input)
fus_fam_simulated<-bootMer(lm_model,function(x) predict(x,data.frame(sang_prop=0),re.form=NA),nsim=simulation_number)$t
pred_values_fusca <- exp(fus_fam_simulated)


summary(pred_values_sang-pred_values_fusca)

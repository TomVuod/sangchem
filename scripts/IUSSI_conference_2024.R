library(sangchem)
library(dplyr)
library(ggplot2)
library(mixOmics)
library(lmerTest)
data("development_data")
data("mass_spectra_data")
dev_samples <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor", "aggression target", "queenless colony")) %>%
  dplyr::select(chromatogram_ID, species, colony, sang_prop, callow, census_date)
normalized_abundances <- peak_proportions_table(mass_spectra_data[mass_spectra_data$chromatogram_ID %in% dev_samples$chromatogram_ID,])
species_indices <- calculate_SII(normalized_abundances)
dev_samples <- left_join(dev_samples, species_indices)



model_input <- dev_samples[dev_samples$species=="F. sanguinea"&dev_samples$callow==0,]
plot_data <- dev_samples[dev_samples$callow==0,]

lm_model <- lmer(predicted_species~sang_prop+(1|colony)+(1|colony:census_date), data=model_input)
x <- seq(0,1,by=0.001)
model_prediction <- data.frame(sang_prop = x, predicted_species=predict(lm_model, data.frame(sang_prop=x), re.form=NA),
                               species = "F. sanguinea")
model_input <- dev_samples[dev_samples$species=="F. fusca",]
lm_model <- lmer(predicted_species~sang_prop+(1|colony)+(1|colony:census_date), data=model_input)



model_prediction <- rbind(model_prediction, data.frame(sang_prop = x, predicted_species=predict(lm_model, data.frame(sang_prop=x), re.form=NA),
                               species = "F. fusca"))

ggplot(plot_data, aes(x=sang_prop, y=predicted_species))+
  geom_point(aes(color=species))+
  xlab("Proportion of the F. sanguinea ants in a colony")+
  ylab("Species identity index")+
  geom_line(aes(color = species),data = model_prediction)







dev_samples <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor", "aggression target", "queenless colony"), callow==0) %>%
  dplyr::select(chromatogram_ID, species, colony, sang_prop, head_width, mass, census_date) %>%
  mutate(corrected_mass=mass/head_width^2*10^6)
model_input <- dev_samples[dev_samples$species=="F. sanguinea",]
# remove outlier
model_input <- model_input[-which.max(model_input$corrected_mass),]
lm_model <- lmer(I(log(corrected_mass)) ~ sang_prop  + (1|colony) +(1 | colony:census_date), data=model_input)
x <- seq(0,1,by=0.001)
model_prediction <- data.frame(sang_prop = x, predicted_mass=exp(predict(lm_model, data.frame(sang_prop=x), re.form=NA)),
                               species = "F. sanguinea")


plot_data <- dev_samples
# remove outlier
plot_data <- plot_data[-which.max(plot_data$corrected_mass),]


model_input <- dev_samples[dev_samples$species=="F. fusca",]

lm_model <- lmer(I(log(corrected_mass)) ~ sang_prop  + (1|colony), data=model_input)
model_prediction <- rbind(model_prediction, data.frame(sang_prop = x, predicted_mass=exp(predict(lm_model, data.frame(sang_prop=x), re.form=NA)),
                               species = "F. fusca"))

ggplot(plot_data, aes(x=sang_prop, y=corrected_mass, color=species, fill=species))+
  scale_color_manual(values=c("black", "red"))+
  geom_point()+
  xlab("Proportion of the F. sanguinea ants in a colony")+
  ylab("Corrected CHC maass")+
  geom_line(aes(y= predicted_mass),data = model_prediction, lwd=.9)




temp<-filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor", "aggression target", "queenless colony")) %>%
  dplyr::select(chromatogram_ID, species, colony, sang_prop, callow, census_date, mass, head_width)
temp<-filter(temp,species=="F. sanguinea")
temp$stand_mass<-(temp$mass/(temp$head_width/1300)^2)
temp$callow[temp$callow==1]<-"callow"
temp$callow[temp$callow!="callow"]<-"mature"



plot_3<-ggplot(temp,aes(y = stand_mass,x=callow)) +
  geom_boxplot()+
  geom_jitter(width=0.2,size=2,col="#0c59d9")+
  scale_y_continuous(trans="log2",breaks=2^(0:5))+
  labs(x="Age",y=expression(paste("Cucitular hydrocarbons mass per individual [",mu,"g]")), title="Cuticular hydrocarbons amount\n in mature and callow F. sanguinea ants")+
  theme(plot.title = element_text(hjust = 0.5),axis.title = element_text(size=14),
        axis.text = element_text(size=12))

# ggsave(file="callow_mature_CHCs_mass.pdf", plot=plot_3, width=7, height=5.8)

data("distances_to_free_living_fusca")
data("separation_experiment_data")

treatments <- unique(separation_experiment_data$treatment)[1:4]

plot_data <- data.frame()
for(i in 1:4){
  plot_data <- rbind(plot_data, data.frame(distance=distances_to_free_living_fusca[[1]][[i]]$rand_distances, treatment=treatments[i]))
}


plot_data_2 <- data.frame()
for(i in 1:4){
  plot_data_2 <- rbind(plot_data_2, data.frame(distance=distances_to_free_living_fusca[[1]][[i]]$true_distance, treatment=treatments[i]))
}


plot_data$treatment <- factor(plot_data$treatment, levels = treatments)


ggplot(plot_data, aes(y=distance, x=treatment))+
  ggdist::stat_halfeye(justification = -0.2, fill = "#cca22c")+
  geom_boxplot(width=0.3, color = "#cca22c", lwd=0.9)+
  xlab("Adult workers age")+
  ylab("Chemical distance to the free-living F. fusca")+
  scale_y_continuous(limits = c(0.28,0.82), breaks = seq(0.3,0.7, 0.1))+
  theme(axis.title=element_text(size = 14),
        axis.text = element_text(size=11),
        panel.background = element_blank(),
          panel.grid.major.y =  element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey"))+
  geom_boxplot(aes(y=distance, x=treatment), data = plot_data_2, color = "red", fill="red")






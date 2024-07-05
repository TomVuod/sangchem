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
  geom_point(aes(color=species), size = 2)+
  xlab("Proportion of the F. sanguinea ants in a colony")+
  ylab("Species identity index")+
  geom_line(aes(color = species),data = model_prediction)+
  scale_color_manual(values = c("black", "red"))+
  theme(axis.title=element_text(size = 16),
        axis.text = element_text(size=13),
        panel.background = element_blank(),
        panel.grid.major =  element_line(linewidth = 0.5, linetype = 'solid',
                                         colour = "grey"),
        legend.text=element_text(size=14))







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

data("CHC_mass_prediction")

plot_data_2 <- data.frame(corrected_mass = CHC_mass_prediction$pred_values_sang, species = "F. sanguinea",
                          sang_prop = 1)
plot_data_2 <- rbind(plot_data_2, data.frame(corrected_mass = CHC_mass_prediction$pred_values_fusca, species = "F. fusca",
                                             sang_prop = 0))




ggplot(plot_data, aes(x=sang_prop, y=corrected_mass, color=species, fill=species))+
  scale_color_manual(values=c("black", "red"))+
  geom_point()+
  xlab("Proportion of the F. sanguinea ants in a colony")+
  ylab("Corrected CHC maass")+
  geom_line(aes(y= predicted_mass),data = model_prediction, lwd=.9)+
  geom_boxplot(data=plot_data_2, aes(x=sang_prop), width = 0.1, fill="#00000000")+
  scale_x_continuous(limits = c(-0.1, 1.1), breaks = seq(0,1,0.1))+
  scale_y_continuous(limits = c(0, 11.5))+
  theme(axis.title=element_text(size = 16),
        axis.text = element_text(size=13),
        panel.background = element_blank(),
        panel.grid.major =  element_line(linewidth = 0.5, linetype = 'solid',
                                         colour = "grey"),
        legend.text=element_text(size=14))

ggplot(plot_data, aes(x=sang_prop, y=corrected_mass, color=species, fill=species))+
  scale_color_manual(values=c("black", "red"))+
  geom_point()+
  xlab("Proportion of the F. sanguinea ants in a colony")+
  ylab("Corrected CHC maass")+
  geom_line(aes(y= predicted_mass),data = model_prediction, lwd=.9)+
  scale_x_continuous(limits = c(-0.1, 1.1), breaks = seq(0,1,0.1))+
  scale_y_continuous(limits = c(0, 11.5))+
  theme(axis.title=element_text(size = 16),
        axis.text = element_text(size=13),
        panel.background = element_blank(),
        panel.grid.major =  element_line(linewidth = 0.5, linetype = 'solid',
                                         colour = "grey"),
        legend.text=element_text(size=14))




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
  theme(axis.title=element_text(size = 16),
        axis.text = element_text(size=13),
        panel.background = element_blank(),
        panel.grid.major.y =  element_line(linewidth = 1, linetype = 'solid',
                                           colour = "grey"))+
  geom_boxplot(aes(y=distance, x=treatment), data = plot_data_2, color = "red", fill="red")







peak_color<-function(peak_ID){
  if (is.na(match(peak_ID,as.numeric(names(species_correlations)))))
    return(rgb(1,1,1))
  corr<-species_correlations[as.character(peak_ID)]
  if (corr>0){
    
    return(rgb(0.5*corr+0.5,-0.5*corr+0.5,-0.5*corr+0.5))
  }
  rgb(0.5*corr+0.5,-0.5*corr+0.5,0.5*corr+0.5)
}


plot_bars<-function(data,descr="",czynnik="density_peak",peak_set=1:110,data2=NULL,proportions=FALSE){
  if (nrow(data)==0) stop ("Empty dataset 1 passed")
  if (!is.null(data2)){
    if(nrow(data2)==0) {
      data2<-NULL 
      warning("Empty dataset 2 passed")
    }
  }
  if (!is.null(data2)&length(descr)!=2) stop("Character vector of the length two should be passed as \"descr\" parameter")
  library(dplyr)
  all_data<-list(data)
  if (!is.null(data2)) all_data<-c(all_data,list(data2))
  if (!proportions){
    for (i in 1:length(all_data)){
      data_table<-all_data[[i]]
      data_table<-data_table[!is.na(data_table[,czynnik]),]
      data_table<-data_table[is.na(match(data_table[,czynnik],c(0,1))),]
      total_abundance<-sum(data_table$abundance)
      data_table$proportion<-data_table$abundance/total_abundance
      data_table<-mutate(data_table,proportion=abundance/sum(abundance))
      data_table<-group_by(data_table,density_peak)#should be changed
      data_table<-summarise(data_table,proportion=sum(proportion),characteristic=characteristic[1])
      all_data[[i]]<-ungroup(data_table)
    }
  }
  data<-all_data[[1]]
  if (!is.null(data2)) data2<-all_data[[2]]
  plot.new()
  par(mar=rep(0.5,4))
  if (!is.null(data2)) plot.window(c(0,111), c(-7, 7))
  else plot.window(c(0,111), c(-1, 7))
  lines(c(0,111),c(0,0))
  text(55,5,descr[1],cex=1)
  for (i in peak_set){
    if (sum(data[,czynnik]==i)!=0){
      #fill_col=list(rgb(0.7,0.7,0.7),"green","red")[[as.numeric(data[data[,czynnik]==i,"characteristic"])+1]]
      #if (length(fill_col)==0) fill_col=rgb(0.7,0.7,0.7)
      #if (is.na(fill_col)) fill_col=rgb(0.7,0.7,0.7)
      fill_col=publication$peak_color(i)
      rect(i-0.3,0,i+0.3,15*data[data[,czynnik]==i,"proportion"],col=fill_col)
      text(i,15*data[data[,czynnik]==i,"proportion"]+0.2+(i%%3)*0.2,as.character(i),cex=0.6,srt=90)
    }
  }
  if (!is.null(data2)){
    text(55,-5,descr[2],cex=1)
    for (i in peak_set){
      if (sum(data2[,czynnik]==i)!=0){
        #fill_col=list(rgb(0.7,0.7,0.7),"green","red")[[as.numeric(data[data[,czynnik]==i,"characteristic"])+1]]
        #if (length(fill_col)==0) fill_col=rgb(0.7,0.7,0.7)
        #if (is.na(fill_col)) fill_col=rgb(0.7,0.7,0.7)
        fill_col=publication$peak_color(i)
        rect(i-0.3,0,i+0.3,-15*data2[data2[,czynnik]==i,"proportion"],col=fill_col)
        text(i,-15*data2[data2[,czynnik]==i,"proportion"]-0.2-(i%%3)*0.2,as.character(i),cex=0.6,srt=90)
      }
    }
  }
}



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
# replace zeros with a small number before transformation
peak_prop[peak_prop==0] <- 10^-16
# apply central log ratio transformation
peak_transformed <- t(apply(peak_prop, 1, clr_transformation))

PCA_species_discr <- prcomp(peak_transformed, scale. = TRUE)

#define groups to be compared
Y <- rownames(PCA_species_discr$x) %in% sanguinea_samples

set.seed(9230)
# run the main model
species_discrimination_PC.plsda <- plsda(PCA_species_discr$x, Y, ncomp=10, scale=FALSE)
spec_discr_PC.plsda.perf<-perf(species_discrimination_PC.plsda, nrepeat = 10,folds=5,
                               progressBar = FALSE,validation = "Mfold",
                               auc=TRUE, scale=FALSE)






















suppressMessages({
  library(sangchem)
  library(dplyr)
  library(mixOmics)
})

# load data
data("development_data")
data("mass_spectra_data")

DA_data <- dplyr::filter(development_data, caste=="worker", !remarks %in% c("aggression actor", "aggression target", "queenless colony")) %>%
  dplyr::select(colony, species, callow, census_date, chromatogram_ID, sang_prop) %>%
  split(list(.$colony, .$census_date), drop=TRUE) %>%
  lapply(function(x){if(sum(x$callow)==0) return(NULL); x}) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  purrr::map2(seq_along(.), function(x,y){x$group <- y; x}) %>%
  {function(x) do.call(rbind, x)}()

Y <- (DA_data$callow==1) + (DA_data$species=="F. sanguinea")

peak_prop <- peak_proportions_table(mass_spectra_data[mass_spectra_data$chromatogram_ID %in% DA_data$chromatogram_ID,])
peak_prop <- peak_prop[match(DA_data$chromatogram_ID, rownames(peak_prop)),]
# replace zeros with a small number before transformation
peak_prop[peak_prop==0] <- 10^-16
# apply central log ratio transformation
peak_transformed <- t(apply(peak_prop, 1, clr_transformation))
PCA_callow_discr <- prcomp(peak_transformed, scale. = TRUE)
PCA_callow_discr$x <- withinVariation(PCA_callow_discr$x, data.frame(DA_data$group))

set.seed(1342)
discr_analysis_callow <- splsda(PCA_callow_discr$x, Y=Y ,multilevel=DA_data$group, 
                                ncomp=7, scale = FALSE)
perf.discr_analysis_callow <- perf(discr_analysis_callow, 
                                   validation = "Mfold", 
                                   folds = 4, 
                                   progressBar = FALSE, 
                                   auc = TRUE, 
                                   nrepeat = 30,
                                   scale = FALSE)
plot(perf.discr_analysis_callow, col = color.mixo(1:3), sd = TRUE, legend.position = "horizontal")
selected_n <- 5
discr_analysis_callows_tuned <- tune.splsda(PCA_callow_discr$x,
                                            Y=Y,
                                            nrepeat = 30,
                                            folds = 4,
                                            ncomp = perf.discr_analysis_callow$choice.ncomp[1,1],
                                            test.keepX = seq(1:120), 
                                            dist = 'max.dist',
                                            validation = 'Mfold',
                                            measure = "BER",
                                            progressBar = FALSE,
                                            scale = FALSE) 

discr_analysis_callows_res <- splsda(PCA_callow_discr$x,Y=Y,
                                     multilevel= DA_data$group,
                                     keepX=discr_analysis_callows_tuned$choice.keepX,
                                     ncomp=selected_n,
                                     scale = FALSE)
plotIndiv(discr_analysis_callows_res,
          legend=TRUE, ellipse=TRUE, comp=1:2, col = c("black", "red", "orange"), title = "Discriminant analysis")



data("discrimination_model_validation")


df <- data.frame(AUC = discrimination_model_validation$AUC_observed, 
                 model = "true sample labels")

df <- rbind(df, data.frame(AUC = discrimination_model_validation$AUC_null, 
                 model = "permuted sample labels"))
ggplot(df, aes(y = AUC, x = model))+
  geom_boxplot(fill="#ffb800", alpha = 0.2)+
  ylab("Area under curve (AUC) of receiver operating characteristic (ROC)")+
  xlab("Model")+  
  theme(axis.title=element_text(size = 14),
                       axis.text = element_text(size=11),
                       panel.background = element_blank(),
                       panel.grid.major.y =  element_line(linewidth = 1, linetype = 'solid',
                                                          colour = "grey"))
  



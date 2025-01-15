library(ggrepel)
library(ggplot2)
library(dplyr)
library(sangchem)
library(forcats)
library(ggpubr)
library(rstatix)
library(stringr)
library(HDInterval)



load_globals()
data("development_data")
paired_observations <- function(group_1, group_2, values_column, group_var="callow", peak_subset=NULL){
  df <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor",
                                                                                     "aggression target",
                                                                                     "queenless colony")) %>%
    select(colony, census_date, chromatogram_ID,
           callow, mass, head_width, species, sang_prop)
  if(values_column=="mass")
    df <- mutate(df, mass=mass/(head_width/1300)^2)
  else
    df$proportion <- subsample_proportion(peak_subset, df$chromatogram_ID)

  df1 <- filter(df, species==!!group_1$species, callow==group_1$callow) %>%
    select(colony, census_date, chromatogram_ID, callow, !!group_var, sang_prop, !!values_column)

  # duplicated samples should converge, thus they should have the same position on x axis
  # each level of the factor correspond to the same colony/sampling date combination
  df1$point_x_pos <- fct_cross(df1$colony, as.character(df1$census_date))
  levels(df1$point_x_pos) <- sample(seq(-0.15,0.15,length.out = length(levels(df1$point_x_pos))))
  df1$point_x_pos <- as.numeric(as.character(df1$point_x_pos))

  filter(df, species==!!group_2$species, callow==group_2$callow) %>%
    select(colony, census_date, chromatogram_ID, callow, !!group_var, sang_prop, !!values_column) -> df2
  # pair with callow ants form the same colony and sampled at the same time
  df <- inner_join(df1, df2, relationship="many-to-many", by = c("colony", "census_date", "sang_prop"))
  # add pair ID
  df$ID <- as.numeric(fct_cross(as.character(df$chromatogram_ID.x), as.character(df$chromatogram_ID.y)))
  df1 <- unique(df[,grep(".*\\.y$", colnames(df), invert = TRUE)])
  colnames(df1) <- gsub("\\.x", "",colnames(df1))
  df2 <- unique(df[,grep(".*\\.x$", colnames(df), invert = TRUE)])
  colnames(df2) <- gsub("\\.y", "",colnames(df2))
  df<-rbind(df1, df2[,colnames(df1)])
  df[,group_var] <- df[,group_var]==group_1[[group_var]]
  return(df)
}

df <- paired_observations(list(species = "F. sanguinea", callow = 0),
                          group_2=list(species = "F. fusca", callow = 0),
                          peak_subset = peaks_sang,
                          values_column="proportion", group_var="species"
)

df$species <- factor(df$species, levels = c(TRUE, FALSE))
levels(df$species) <- c("F. sanguinea", "F. fusca")
p1 <- ggplot(aes(x=species, y=proportion , group=species), data=df)+
  geom_boxplot(col="#1f4a24",
               fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(fill=sang_prop), shape=21,cex=3.4,color="black", position = position_jitter(width=0.1))+
  #geom_line(aes(x=callow+point_x_pos, y=proportion, group=ID), col="#444444")+
  #scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("Species")+
  ylab("Proportion of the F. sanguinea markers\nin the total CHC mass")+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=18)) +
  theme(axis.title.y=element_text(vjust = 0.5))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=13),
        legend.position = "none")+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  geom_signif(
    y_position = max(df[["proportion"]])*1.1, xmin = 1, xmax = 2,
    annotation = "NS", tip_length = 0.02)


df <- paired_observations(list(species = "F. sanguinea", callow = 0),
                          group_2=list(species = "F. fusca", callow = 0),
                          peak_subset = peaks_fusca,
                          values_column="proportion", group_var="species"
)

df$species <- factor(df$species, levels = c(TRUE, FALSE))
levels(df$species) <- c("F. sanguinea", "F. fusca")
p2 <- ggplot(aes(x=species, y=proportion , group=species), data=df)+
  geom_boxplot(col="#1f4a24",
               fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(fill=sang_prop), shape=21,cex=3.4,color="black", position = position_jitter(width=0.1))+
  #geom_line(aes(x=callow+point_x_pos, y=proportion, group=ID), col="#444444")+
  #scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("Species")+
  ylab("Proportion of the F. fusca markers\nin the total CHC mass")+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=18)) +
  theme(axis.title.y=element_text(vjust = 0.5))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=13), legend.position = "none")+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  geom_signif(
    y_position = max(df[["proportion"]])*1.1, xmin = 1, xmax = 2,
    annotation = "NS", tip_length = 0.02)




df <- paired_observations(list(species = "F. sanguinea", callow = 0),
                          group_2=list(species = "F. sanguinea", callow = 1),
                          values_column="mass", group_var="callow")
df$callow <- as.numeric(df$callow)

p3 <- ggplot(aes(x=callow, y=mass, group=callow), data=df)+
  geom_boxplot(col="#1f4a24",
               fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(x=callow+point_x_pos, y=mass, fill=sang_prop), shape=21, cex=3.4, color="black")+
  geom_line(aes(x=callow+point_x_pos, y=mass, group=ID), col="#444444")+
  scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_x_continuous(breaks = c(0,1), labels = c("callow", "mature"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea - age")+
  ylab(expression(paste("Mass per individual [",mu,"g]")))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=18))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=13))+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"),
        legend.position = "none")+
  geom_signif(
    y_position = log2(max(df[["mass"]])*1.1), xmin = 0, xmax = 1,
    annotation = 0.0000305, tip_length = 0.02
  )






n_alkanes_ID <- unique(mass_spectra_data$peak_ID[grep("^C[0-9]{2}$",mass_spectra_data$identification)])

df <- paired_observations(list(species = "F. sanguinea", callow = 0),
                          group_2=list(species = "F. sanguinea", callow = 1),
                          peak_subset = peaks_sang,
                          values_column="proportion", group_var="callow"
)

df$mature<- as.numeric(df$callow)
p6 <- ggplot(aes(x=as.numeric(callow), y=proportion , group=mature), data=df)+
  geom_boxplot(col="#1f4a24",
               fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(x=as.numeric(callow)+point_x_pos, y=proportion, fill=sang_prop),shape=21,cex=3.4,color="black")+
  #geom_line(aes(x=callow+point_x_pos, y=proportion, group=ID), col="#444444")+
  #scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_x_continuous(breaks = c(0,1), labels = c("callow", "mature"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea - age")+
  ylab("Proportion of the n-alkanes in the total CHC mass")+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=18)) +
  theme(axis.title.y=element_text(vjust = 0.5))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=13),
        legend.position = "none")+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  geom_signif(
    y_position = max(df[["proportion"]])*1.1, xmin = 0, xmax = 1,
    annotation = 0.000153, tip_length = 0.02)


df <- paired_observations(list(species = "F. sanguinea", callow = 0),
                          group_2=list(species = "F. sanguinea", callow = 1),
                          peak_subset = n_alkanes_ID,
                          values_column="proportion", group_var="callow"
)

df$mature<- as.numeric(df$callow)
p7 <- ggplot(aes(x=as.numeric(callow), y=proportion , group=mature), data=df)+
  geom_boxplot(col="#1f4a24",
               fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(x=as.numeric(callow)+point_x_pos, y=proportion, fill=sang_prop),shape=21,cex=3.4,color="black")+
  #geom_line(aes(x=callow+point_x_pos, y=proportion, group=ID), col="#444444")+
  #scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_x_continuous(breaks = c(0,1), labels = c("callow", "mature"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea - age")+
  ylab("Proprotion of F. sanguinea markers in the total CHC mass")+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=18)) +
  theme(axis.title.y=element_text(vjust = 0.5))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=13),
        legend.position = c(0.3, 0.7))+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  geom_signif(
    y_position = max(df[["proportion"]])*1.1, xmin = 0, xmax = 1,
    annotation = 0.000153, tip_length = 0.02)


dev.new()
pdf("Fig1.pdf", width = 19, height = 7.5)
ggarrange(p1, p2, p3, p6, p7, ncol=5, labels=c("A", "B", "C", "D", "E"), widths = c(4.5, 4.5, 4, 4,4))
dev.off()



# Figure 2

data("development_data")
dev_samples <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor", "aggression target", "queenless colony")) %>%
  dplyr::select(chromatogram_ID, species, colony, sang_prop, callow, census_date)
normalized_abundances <- peak_proportions_table(mass_spectra_data[mass_spectra_data$chromatogram_ID %in% dev_samples$chromatogram_ID,])
species_indices <- calculate_SII(normalized_abundances)
dev_samples <- left_join(dev_samples, species_indices)
model_input <- dev_samples[dev_samples$species=="F. sanguinea"&dev_samples$callow==1,]
model_input$SII_difference <- model_input$predicted_species - find_nestmates_SII(model_input$chromatogram_ID, dev_samples, heterospecific = FALSE)
# no random effect since it captures no variance
lm_model <- lm(SII_difference~sang_prop, data=model_input)
x <- seq(0,1,by=0.001)
model_prediction <- data.frame(sang_prop = x, SII_difference=predict(lm_model, data.frame(sang_prop=x), re.form=NA))

dev.new()
pdf("Fig2.pdf", width = 7, height = 4)
ggplot(model_input, aes(x=sang_prop, y=SII_difference))+
  geom_point(size = 2)+
  xlab("Proportion of the F. sanguinea ants in a colony")+
  ylab("Difference between Species Identity Indices")+
  geom_line(data = model_prediction, lwd=1)+
  scale_color_manual(values = c("black", "red"))+
  scale_y_continuous(limits = c(-0.44, 0.25))+
  theme(axis.title=element_text(size = 16),
        axis.text = element_text(size=13),
        panel.background = element_blank(),
        panel.grid.major =  element_line(linewidth = 0.5, linetype = 'solid',
                                         colour = "grey"),
        legend.text=element_text(size=14))+
  geom_hline(aes(yintercept = 0), lwd=1.1, linetype = "dashed", color = "#e15a00")
dev.off()




# Figure 3

compare_age_groups <- function(x){
  x$peak_number <- unlist(lapply(x$chromatogram_ID, function(x) nrow(mass_spectra_data[mass_spectra_data$chromatogram_ID==x & mass_spectra_data$peak_ID > 1,])))
  x
}

differential_profile <- function(x){
  callow_IDs <- x[x$callow==1, "chromatogram_ID"]
  mature_IDs <- x[x$callow==0, "chromatogram_ID"]
  if(length(callow_IDs)==0 | length(mature_IDs)==0) return(NULL)
  print(head(mass_spectra_data))
  average_profile_callow <- averaged_profile(callow_IDs, mass_spectra_data)
  average_profile_mature <- averaged_profile(mature_IDs, mass_spectra_data)
  diff_averaged_profile <- merge(average_profile_callow, average_profile_mature, all=TRUE, by="peak_ID")
  diff_averaged_profile[,"relative_abundance.x"][is.na(diff_averaged_profile[,"relative_abundance.x"])] <- 0
  diff_averaged_profile[,"relative_abundance.y"][is.na(diff_averaged_profile[,"relative_abundance.y"])] <- 0
  print("end")
  data.frame(peak_ID = diff_averaged_profile[,"peak_ID"], relative_abundance = diff_averaged_profile[,"relative_abundance.x"] - diff_averaged_profile[,"relative_abundance.y"])
}


suppressMessages({
  library(sangchem)
  library(dplyr)
  library(mixOmics)
})

# load data
data("development_data")
data("mass_spectra_data")

DA_data <- dplyr::filter(development_data, caste=="worker", !remarks %in% c("aggression actor", "aggression target", "queenless colony"), species=="F. sanguinea") %>%
  dplyr::select(colony, species, callow, census_date, chromatogram_ID, sang_prop) %>%
  split(list(.$colony, .$census_date), drop=TRUE) %>%
  lapply(function(x){if(sum(x$callow)==0) return(NULL); x}) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  lapply(differential_profile) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  purrr::reduce(function(x,y) merge(x,y,all=TRUE,by="peak_ID"))

DA_data <- as.data.frame(lapply(as.list(DA_data), function(x) {x[is.na(x)] <- 0; x}))

colnames(DA_data)[2:length(colnames(DA_data))] <- paste0("col_", 1:(length(colnames(DA_data))-1))
diff_data <- tidyr::gather(DA_data, key = "sample", value = "relative_abundance", -peak_ID)

age_profile_sample <- function(x){
  callow_IDs <- x[x$callow==1, "chromatogram_ID"]
  mature_IDs <- x[x$callow==0, "chromatogram_ID"]
  if(length(callow_IDs)==0 | length(mature_IDs)==0) return(NULL)
  average_profile_callow <- averaged_profile(callow_IDs, mass_spectra_data)
  average_profile_mature <- averaged_profile(mature_IDs, mass_spectra_data)
  average_profile_callow$callow <- 1
  average_profile_mature$callow <- 0
  rbind(average_profile_callow,  average_profile_mature)
}



age_profiles <- dplyr::filter(development_data, caste=="worker", !remarks %in% c("aggression actor", "aggression target", "queenless colony"), species=="F. sanguinea") %>%
  dplyr::select(colony, species, callow, census_date, chromatogram_ID, sang_prop) %>%
  split(list(.$colony, .$census_date), drop=TRUE) %>%
  lapply(function(x){if(sum(x$callow)==0) return(NULL); x}) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  lapply(age_profile_sample) %>%
  {function(x) x[!unlist((lapply(x, is.null)))]}() %>%
  purrr::reduce(rbind) %>%
  split(.$callow, drop=TRUE) %>%
  lapply(mean_profile)


all_peaks <- unique(c(age_profiles[[1]]$peak_ID, age_profiles[[2]]$peak_ID))
missing_peaks <- setdiff(all_peaks, age_profiles[[1]]$peak_ID)
if(length(missing_peaks)>0)
  age_profiles[[1]] <- rbind(age_profiles[[1]], data.frame(peak_ID = missing_peaks, relative_abundance = 0))
age_profiles[[1]]$age = "mature"
missing_peaks <- setdiff(all_peaks, age_profiles[[2]]$peak_ID)
if(length(missing_peaks)>0)
  age_profiles[[2]] <- rbind(age_profiles[[2]], data.frame(peak_ID = missing_peaks, relative_abundance = 0))
age_profiles[[2]]$age = "callow"
#plot_data <- rbind(age_profiles[[1]], age_profiles[[2]])


devtools::load_all()
source(system.file("/scripts/peak_ID_publication.R", package="sangchem"))
age_profiles[[1]]$peak_ID <- match(age_profiles[[1]]$peak_ID,peaks_ordered)
age_profiles[[2]]$peak_ID <- match(age_profiles[[2]]$peak_ID,peaks_ordered)
diff_data$peak_ID <- match(diff_data$peak_ID, peaks_ordered)


library(ggplot2)

dev.new()
pdf("Fig3.pdf", width = 12, height = 6)
ggplot(NULL, aes(x=peak_ID, y=relative_abundance, fill = age, color=age))+
  geom_bar(data = age_profiles[[1]], stat = "identity",linewidth=0.2,
           position = position_dodge(1), width = 0.6)+
  geom_bar(data = age_profiles[[2]], stat = "identity",linewidth=0.2,
           position = position_dodge(1), width = 0.6)+
  scale_color_manual(values=c("transparent", "#00000099"))+
  scale_fill_manual(values=c("#99999966", "transparent"))+
  geom_text(data = age_profiles[[1]], aes(label=peak_ID), size = 2,
            hjust=-0.8, color="black", angle = 90, vjust = 0.3)+
  geom_point(data = diff_data, aes(x=peak_ID, y=relative_abundance), color = "black", fill="black", size = 0.1)+
  theme(panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_rect(color="white", fill = "white"))+
  ylab("Relative abundance")+
  xlab("Peak ID")

dev.off()



# Figure 4

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



dev.new()
pdf("Fig4.pdf", width = 9, height = 5)
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

dev.off()





# Figure 5

data("distances_to_free_living_fusca")

#' @importFrom HDInterval hdi
plot_data <- data.frame()
for(i in 1:length(distances_to_free_living_fusca$all)){
  distances <- distances_to_free_living_fusca$all[[i]]$rand_distance
  dist_density <- density(distances, n=10^3)
  plot_data_part <- data.frame(x = dist_density$x, y = dist_density$y)
  interval <- hdi(distances)
  plot_data_part$treatment <- names(distances_to_free_living_fusca$all)[i]
  plot_data_part$HDI_mask <- FALSE
  plot_data_part$true_distance <- distances_to_free_living_fusca$all[[i]]$true_distance
  plot_data_part$max_y <- max(plot_data_part$y)
  plot_data_part$HDI_mask[plot_data_part$x>interval[1] & plot_data_part$x<interval[2]] <- TRUE
  plot_data <- rbind(plot_data, plot_data_part)
}


dev.new()
pdf("Fig5.pdf", width = 5, height = 7)
ggplot(plot_data, aes(x=x, y=y))+
  facet_wrap(~treatment, labeller = labeller(treatment=c(treatment_1 = "1-3 days", treatment_2 = "7-10 days",
                                                         treatment_3 = "17-20 days", treatment_4 = "35-40 days")),
             ncol=1)+
  geom_segment(aes(x = true_distance, y = 0, xend = true_distance, yend = max_y), lwd=1, lty=2, color="#cc0033")+
  geom_area(aes(x = x, y = y),
            data = subset(plot_data, HDI_mask), alpha = 0.5, fill = "#11bb22")+
  geom_line()+
  ylab("Density")+
  xlab("Chemical distance")+
  theme(panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_rect(color="black", fill = "white"))
dev.off()
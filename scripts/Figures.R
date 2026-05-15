library(ggrepel)
library(ggplot2)
library(sangchem)
library(forcats)
library(ggpubr)
library(rstatix)
library(stringr)
library(HDInterval)
library(mixOmics)
library(dplyr)
library(lme4)



load_globals()
data("development_data")

data("body_surface_area")
CHC_normalizer <- get_body_area_normalizer(body_surface_area, development_data)

paired_observations <- function(group_1, group_2, values_column, group_var="callow", peak_subset=NULL){
  df <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor",
                                                                                     "aggression target",
                                                                                     "queenless colony")) %>%
    select(colony, census_date, chromatogram_ID,
           callow, mass, head_width, species, sang_prop)
  if(values_column=="mass")
    df <- mutate(df, mass=mass*CHC_normalizer(species, head_width))
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
  geom_boxplot(col="#633D0F",
               fill="#C2C7BD66", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(fill=sang_prop), shape=21,cex=2.7,color="black", position = position_jitter(width=0.14))+
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
        legend.position = c(0.77, 0.81),
        legend.background = element_rect(fill = NA))+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  geom_signif(
    y_position = max(df[["proportion"]])*1.1, xmin = 1, xmax = 2,
    annotation = "6.1e-5", tip_length = 0.02)


df <- paired_observations(list(species = "F. sanguinea", callow = 0),
                          group_2=list(species = "F. fusca", callow = 0),
                          peak_subset = peaks_fusca,
                          values_column="proportion", group_var="species"
)

df$species <- factor(df$species, levels = c(TRUE, FALSE))
levels(df$species) <- c("F. sanguinea", "F. fusca")
p2 <- ggplot(aes(x=species, y=proportion , group=species), data=df)+
  geom_boxplot(col="#633D0F",
               fill="#C2C7BD66", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(fill=sang_prop), shape=21,cex=2.7,color="black", position = position_jitter(width=0.14))+
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
    annotation = "6.1e-5", tip_length = 0.02)




df <- paired_observations(list(species = "F. sanguinea", callow = 0),
                          group_2=list(species = "F. fusca", callow = 0),
                          values_column="mass", group_var="species"
)

df$species <- factor(df$species, levels = c(TRUE, FALSE))
levels(df$species) <- c("F. sanguinea", "F. fusca")

p3 <- ggplot(aes(x=species, y=mass , group=species), data=df)+
  geom_boxplot(col="#633D0F",
               fill="#C2C7BD66", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(fill=sang_prop), shape=21,cex=2.7,color="black", position = position_jitter(width=0.14))+
  #geom_line(aes(x=callow+point_x_pos, y=proportion, group=ID), col="#444444")+
  #scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("Species")+
  ylab(expression(paste("Normalized CHC mass [",mu,"g]")))+
  scale_y_continuous(trans="log2", breaks=2^(0:5))+
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
    y_position = 5, xmin = 1, xmax = 2,
    annotation = "0.00043", tip_length = 0.02)





df <- paired_observations(list(species = "F. sanguinea", callow = 0),
                          group_2=list(species = "F. sanguinea", callow = 1),
                          values_column="mass", group_var="callow")
df$callow <- as.numeric(df$callow)

p4 <- ggplot(aes(x=callow, y=mass, group=callow), data=df)+
  geom_boxplot(col="#633D0F",
               fill="#C2C7BD66", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(fill=sang_prop), shape=21,cex=2.7,color="black", position = position_jitter(width=0.14))+
  geom_line(aes(x=callow+point_x_pos, y=mass, group=ID), col="#444444")+
  scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_x_continuous(breaks = c(0,1), labels = c("callow", "mature"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea - age")+
  ylab(expression(paste("Normalized CHC mass [",mu,"g]")))+
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
                          peak_subset = n_alkanes_ID,
                          values_column="proportion", group_var="callow"
)

df$mature<- as.numeric(df$callow)
p6 <- ggplot(aes(x=as.numeric(callow), y=proportion , group=mature), data=df)+
  geom_boxplot(col="#633D0F",
               fill="#C2C7BD66", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(fill=sang_prop), shape=21,cex=2.7,color="black", position = position_jitter(width=0.14))+
  #geom_line(aes(x=callow+point_x_pos, y=proportion, group=ID), col="#444444")+
  #scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_x_continuous(breaks = c(0,1), labels = c("callow", "mature"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea - age")+
  ylab("Proportion of n-alkanes in the total CHC mass")+
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
    annotation = 0.0000305, tip_length = 0.02)




dev.new()
pdf("./inst/figures/Fig2.pdf", width = 19, height = 7.5)
#dev.new()
ggarrange(p1, p2, p4, p3, p6, ncol=5, labels=LETTERS[1:5], widths = c(4.5, 4.5, 4, 4, 4))
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
  average_profile_callow <- averaged_profile(callow_IDs, mass_spectra_data)
  average_profile_mature <- averaged_profile(mature_IDs, mass_spectra_data)
  diff_averaged_profile <- merge(average_profile_callow, average_profile_mature, all=TRUE, by="peak_ID")
  diff_averaged_profile[,"relative_abundance.x"][is.na(diff_averaged_profile[,"relative_abundance.x"])] <- 0
  diff_averaged_profile[,"relative_abundance.y"][is.na(diff_averaged_profile[,"relative_abundance.y"])] <- 0
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

# new peak IDs
batch7 <- mass_spectra_data[mass_spectra_data$batch==7 & mass_spectra_data$peak_ID!=0,]
peaks <- unique(batch7$peak_ID)
mean_ret_times <- unlist(lapply(peaks, function(x) mean(batch7[batch7$peak_ID==x, "retention_time"])))
peaks_ordered <- peaks[order(mean_ret_times)]

peak_ID_df <- data.frame(peak_ID = 1:length(peaks_ordered), old_ID = peaks_ordered)
age_profiles[[1]]$peak_ID <-  peak_ID_df$peak_ID[match(age_profiles[[1]]$peak_ID, peak_ID_df$old_ID)]
age_profiles[[2]]$peak_ID <- peak_ID_df$peak_ID[match(age_profiles[[2]]$peak_ID, peak_ID_df$old_ID)]
diff_data$peak_ID <- peak_ID_df$peak_ID[match(diff_data$peak_ID, peak_ID_df$old_ID)]


library(ggplot2)

dev.new()
pdf("./inst/figures/Fig3.pdf", width = 8, height = 4)
ggplot(NULL, aes(x=peak_ID, y=relative_abundance, fill = age, color=age))+
  geom_bar(data = age_profiles[[1]], stat = "identity",linewidth=0.2,
           position = position_dodge(1), width = 0.6)+
  geom_bar(data = age_profiles[[2]], stat = "identity",linewidth=0.2,
           position = position_dodge(1), width = 0.6)+
  scale_color_manual(values=c("transparent", "#00000099"))+
  scale_fill_manual(values=c("#99999966", "transparent"))+
  geom_text(data = age_profiles[[1]], aes(label=peak_ID), size = 2,
            hjust=-2, color="black", angle = 90, vjust = 0.3)+
  geom_point(data = diff_data, aes(x=peak_ID, y=relative_abundance), color = "black", fill="black", size = 0.1)+
  theme(panel.border = element_rect(color = "black", fill=NA),
        panel.background = element_rect(color="white", fill = "white"),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank())+
  ylab("Relative abundance")+
  xlab("Peak ID")

dev.off()



# Figure 4

dev_samples <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor", "aggression target", "queenless colony"), callow==0) %>%
  dplyr::select(chromatogram_ID, species, colony, sang_prop, head_width, mass, census_date) %>%
  mutate(normalized_mass=mass*CHC_normalizer(species, head_width))
model_input <- dev_samples[dev_samples$species=="F. sanguinea",]
# remove outlier
model_input <- model_input[-which.max(model_input$normalized_mass),]
lm_model <- lmer(I(log(normalized_mass)) ~ sang_prop+(1|colony) +(1 | colony:census_date), data=model_input)
x <- seq(0,1,by=0.001)
model_prediction <- data.frame(sang_prop = x, predicted_mass=exp(predict(lm_model, data.frame(sang_prop=x,head_width=1.3), re.form=NA)),
                               species = "F. sanguinea")


plot_data <- dev_samples
# remove outlier
plot_data <- plot_data[-which.max(plot_data$normalized_mass),]

model_input <- dev_samples[dev_samples$species=="F. fusca",]

lm_model <- lmer(I(log(normalized_mass)) ~ sang_prop + (1|colony), data=model_input)
model_prediction <- rbind(model_prediction, data.frame(sang_prop = x, predicted_mass=exp(predict(lm_model, data.frame(sang_prop=x,head_width=1.3), re.form=NA)),
                                                       species = "F. fusca"))

data("CHC_mass_prediction")

plot_data_2 <- data.frame(normalized_mass = CHC_mass_prediction$pred_values_sang, species = "F. sanguinea",
                          sang_prop = 1)
plot_data_2 <- rbind(plot_data_2, data.frame(normalized_mass = CHC_mass_prediction$pred_values_fusca, species = "F. fusca",
                                             sang_prop = 0))



dev.new()
pdf("./inst/figures/Fig4.pdf", width = 9, height = 5)
ggplot(plot_data, aes(x=sang_prop, y=normalized_mass, color=species, fill=species))+
  scale_color_manual(values=c("black", "red"))+
  geom_point(size=2.4, alpha=0.3)+
  xlab("Proportion of the F. sanguinea ants in a colony")+
  ylab(expression(paste("Normalized CHC mass [",mu,"g]")))+
  geom_line(aes(y= predicted_mass),data = model_prediction, lwd=.9)+
  geom_boxplot(data=plot_data_2, aes(x=sang_prop), width = 0.1, fill="#00000000",
               outlier.size = 0.3, outlier.alpha = 1, outlier.shape = 4)+
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

library(dplyr)
library(gridExtra)
library(grid)
library(sangchem)
library(ComplexHeatmap)
plot_list <- list()
separation_periods <- c("1-3 days", "8-10 days", "17-20 days", "35-40 days")
slave_source_colony <- gsub("fus ", "", slave_source_colony)
sslave_source_colony <- gsub("Kam", "K", slave_source_colony)
slave_source_colony <- gsub("Chot", "Ch", slave_source_colony)
names(colony_to_chrID) <- gsub("fus ", "", names(colony_to_chrID))
names(colony_to_chrID) <- gsub("Kam", "K", names(colony_to_chrID))
names(colony_to_chrID) <- gsub("Chot", "Ch", names(colony_to_chrID))


for(n in 1:4){
  treatment_df <- filter(separation_experiment_data, treatment_id ==n)
  treatment_df <- treatment_df[!is.na(colony_to_chrID[slave_source_colony[treatment_df$colony]]),]
  treatment_df <- arrange(treatment_df, colony)
  colonies <- sort(unique(treatment_df$colony))
  source_colonies <- slave_source_colony[colonies]
  source_colonies <- source_colonies[!duplicated(source_colonies)]
  sample_pairs <- expand.grid(treatment_df$chromatogram_ID, colony_to_chrID[source_colonies])
  colnames(sample_pairs) <- c("chr_ID_1", "chr_ID_2")
  dissimilarities <- purrr::pmap_dbl(sample_pairs, pair_distance, MS_data=mass_spectra_data)
  sample_colony <- treatment_df$colony
  sample_colony <- factor(sample_colony, levels = sample_colony[!duplicated(sample_colony)])
  dissimilarity_matrix <- matrix(dissimilarities, ncol = length(source_colonies))
  dissimilarity_matrix <- apply(dissimilarity_matrix, 2, function(x) tapply(x, sample_colony, mean))
  cell_fun = function(j, i, x, y, w, h, fill){
    source_colony <- slave_source_colony[levels(sample_colony)[i]]
    fusca_colony <- source_colonies[j]
    if(source_colony==fusca_colony) {
      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "#2cbf51", lwd=5))
    }
    else{
      grid.rect(x, y, w*0.8, h*0.8, gp = gpar(fill = fill, col = NA))
    }
    grid.text(round(dissimilarity_matrix[i, j],3), x, y, gp = gpar(fontsize = 10))
  }

  if(n<4){
    plot_list[[n]] <- grid.grabExpr(draw(
      ComplexHeatmap::Heatmap(dissimilarity_matrix,
                              col = circlize::colorRamp2(c(0.2, 0.55, 0.9), c("#f74002", "white", "#4d7deb")),
                              name = "Treatment: 1-3 days",
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              column_labels = source_colonies,
                              row_label=levels(sample_colony),
                              cell_fun = cell_fun,
                              row_title=gt_render("<span style='font-size:12pt'>*F. sanguinea* colony ID</span>"),
                              column_title=gt_render(paste0("<span style='font-size:18pt'>",separation_periods[n],"</span><br>",
                                                            "<span style='font-size:12pt'>*F. fusca* colony ID</span>"),
                                                     r = unit(2, "pt"),
                                                     padding = unit(c(2, 2, 2, 2), "pt")),
                              row_names_side = "left",
                              column_names_side = "top",
                              show_heatmap_legend = FALSE,
                              column_names_gp = gpar(fontsize = 9),
                              row_names_gp = gpar(fontsize = 9))

    ))

  }
  else{
    plot_list[[n]] <- grid.grabExpr(draw(
      ComplexHeatmap::Heatmap(dissimilarity_matrix,
                              col = circlize::colorRamp2(c(0.2, 0.55, 0.9), c("#f74002", "white", "#4d7deb")),
                              name = "Treatment: 1-3 days",
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              column_labels = source_colonies,
                              row_label=levels(sample_colony),
                              cell_fun = cell_fun,
                              row_title=gt_render("<span style='font-size:12pt'>*F. sanguinea* colony ID</span>"),
                              column_title=gt_render(paste0("<span style='font-size:18pt'>",separation_periods[n],"</span><br>", "<span style='font-size:12pt'>*F. fusca* colony ID</span>"),  r = unit(2, "pt"),  padding = unit(c(2, 2, 2, 2), "pt")),
                              row_names_side = "left",
                              column_names_side = "top",
                              heatmap_legend_param = list(title = "Chemical\ndistance"),
                              column_names_gp = gpar(fontsize = 9),
                              row_names_gp = gpar(fontsize = 9))
    ))

  }
}

dev.new()
pdf("./inst/figures/Fig5.pdf", width = 8, height = 8)
gridExtra::grid.arrange(plot_list[[1]],
                        plot_list[[2]],
                        plot_list[[3]],
                        plot_list[[4]],ncol=2)
dev.off()

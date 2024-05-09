library(ggrepel)
library(ggplot2)
library(dplyr)
library(sangchem)
library(forcats)
library(ggpubr)
library(rstatix)
library(stringr)



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
        legend.text=element_text(size=13))+
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
        legend.text=element_text(size=13))+
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
  scale_x_continuous(breaks = c(0,1), labels = c("Callow ants", "Mature ants"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea workers")+
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


df <- paired_observations(list(species = "F. sanguinea", callow = 0),
                          group_2=list(species = "F. sanguinea", callow = 1),
                          peak_subset = peaks_callow,
                          values_column="mass", group_var="callow"
)

df$mature<- as.numeric(df$callow)

p4 <- ggplot(aes(x=as.numeric(callow), y=mass , group=mature), data=df)+
  geom_boxplot(col="#1f4a24",
               fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(x=as.numeric(callow)+point_x_pos, y=mass, fill=sang_prop),shape=21,cex=3.4,color="black")+
  geom_line(aes(x=callow+point_x_pos, y=mass, group=ID), col="#444444")+
  #scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_x_continuous(breaks = c(0,1), labels = c("Callow ants", "Mature ants"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea workers")+
  ylab(expression(paste("Mass per individual [",mu,"g]")))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=18)) +
  theme(axis.title.y=element_text(vjust = 0.5))+
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  geom_signif(
    y_position = max(df[["mass"]])*1.1, xmin = 0, xmax = 1,
    annotation = 0.000153, tip_length = 0.02)


df <- paired_observations(list(species = "F. sanguinea", callow = 0),
                          group_2=list(species = "F. sanguinea", callow = 1),
                          peak_subset = peaks_callow,
                          values_column="proportion", group_var="callow"
)

df$mature<- as.numeric(df$callow)
p5 <- ggplot(aes(x=as.numeric(callow), y=proportion , group=mature), data=df)+
  geom_boxplot(col="#1f4a24",
               fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(x=as.numeric(callow)+point_x_pos, y=proportion, fill=sang_prop),shape=21,cex=3.4,color="black")+
  geom_line(aes(x=callow+point_x_pos, y=proportion, group=ID), col="#444444")+
  #scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_x_continuous(breaks = c(0,1), labels = c("Callow ants", "Mature ants"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea workers")+
  ylab("Proportion of the callow markers\nin the total CHC mass")+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=18)) +
        theme(axis.title.y=element_text(vjust = 0.5))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=13),
        legend.position = c(0.8, 0.7))+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  geom_signif(
    y_position = max(df[["proportion"]])*1.1, xmin = 0, xmax = 1,
    annotation = 0.000153, tip_length = 0.02)




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
  scale_x_continuous(breaks = c(0,1), labels = c("Callow ants", "Mature ants"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea workers")+
  ylab("Proportion of the n-alkanes in the total CHC mass")+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=18)) +
  theme(axis.title.y=element_text(vjust = 0.5))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=13),
        legend.position = c(0.8, 0.7))+
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
  scale_x_continuous(breaks = c(0,1), labels = c("Callow ants", "Mature ants"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea workers")+
  ylab("Proprotion of F. sanguinea markers in the total CHC mass")+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=18)) +
  theme(axis.title.y=element_text(vjust = 0.5))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=13),
        legend.position = c(0.8, 0.7))+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  geom_signif(
    y_position = max(df[["proportion"]])*1.1, xmin = 0, xmax = 1,
    annotation = 0.000153, tip_length = 0.02)


ggarrange(p1, p2, p3, p4, p5, ncol=5, labels=c("A", "B", "C", "D", "E"), widths = c(4, 4, 4.3, 4,4))



mean_abundances_diff <- function(df){
  chr_IDs_mature <- filter(df, species =="F. sanguinea", callow==0) %>%
    pull(chromatogram_ID)
  chr_IDs_callow <- filter(df, species =="F. sanguinea", callow==1) %>%
    pull(chromatogram_ID)
  if(length(chr_IDs_mature)==0 | length(chr_IDs_callow)==0) return(NULL)
  mean_peaks_mature <- relative_amounts(mass_spectra_data, chr_IDs_mature)
  mean_peaks_callow <- relative_amounts(mass_spectra_data, chr_IDs_callow)
  mean_peaks_merged <- merge(mean_peaks_callow, mean_peaks_mature, by = "peak_ID")
  mean_peaks_merged$relative_abundance.x[is.na(mean_peaks_merged$relative_abundance.x)] <- 0
  mean_peaks_merged$relative_abundance.y[is.na(mean_peaks_merged$relative_abundance.y)] <- 0
  mean_peaks_merged$relative_abundance <- mean_peaks_merged$relative_abundance.x - mean_peaks_merged$relative_abundance.y
  mean_peaks_merged[,c("peak_ID", "relative_abundance")]
}

df <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor",
                                                                                   "aggression target",
                                                                                   "queenless colony")) %>%
  select(colony, census_date, chromatogram_ID,
         callow, mass, head_width, species, sang_prop)
callow_mature_diff <- split(df, list(df$census_date, df$colony), drop=TRUE) %>%
  lapply(mean_abundances_diff)

callow_mature_diff <- callow_mature_diff[!unlist(lapply(callow_mature_diff, is.null))]
callow_mature_diff <- purrr::reduce(callow_mature_diff, merge, by ="peak_ID", all=TRUE)
for(col_ind in grep("relative_abundance", colnames(callow_mature_diff))){
  callow_mature_diff[,col_ind][is.na(callow_mature_diff[,col_ind])] <- 0
}

plot_data <- data.frame(peak_ID = callow_mature_diff$peak_ID,
                        relative_abundance_diff = apply(callow_mature_diff[,grep("relative_abundance", colnames(callow_mature_diff))], 1, mean),
                        sd = apply(callow_mature_diff[,grep("relative_abundance", colnames(callow_mature_diff))], 1, sd))

# line
ggplot(plot_data, aes(x=peak_ID, y=relative_abundance_diff))+
  geom_bar(stat="identity", fill="grey") +
  geom_text(aes(label = peak_ID))
  #geom_linerange( aes(x=peak_ID, ymin=relative_abundance_diff-sd, ymax=relative_abundance_diff+sd), colour="orange", size=1.3)



data("distances_to_free_living_fusca")

#' @importFrom HDInterval hdi
plot_data <- data.frame()
for(i in 1:length(distances_to_free_living_fusca$all)){
  plot_data_part <- data.frame(distance = distances_to_free_living_fusca$all[[i]]$rand_distance)
  interval <- hdi(plot_data_part$distance)
  plot_data_part$treatment <- names(distances_to_free_living_fusca$all)[i]
  plot_data_part$HDI_mask <- FALSE
  plot_data_part$HDI_mask[plot_data_part$distance>interval[1] & plot_data_part$distance<interval[2]] <- TRUE
  plot_data <- rbind(plot_data, plot_data_part)
}

ggplot(plot_data, aes(x=distance))+
  facet_wrap(~treatment)+
  geom_density()+
  geom_ribbon(
    stat = "density",
    aes(
      ymin = 0,
      ymax = {x <- after_stat(density); x[!plot_data[["HDI_mask"]]] <- 0; x}
    )
  )


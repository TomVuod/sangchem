library(ggrepel)
library(ggplot2)
library(dplyr)
library(sangchem)
library(forcats)
library(ggpubr)
library(rstatix)



load_globals()
data("development_data")
paired_observations <- function(group_1, group_2, values_column, group_var="callow", peak_subset=NULL){
  df <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor",
                                                                                     "aggression target",
                                                                                     "queenless colony"),
               species=="F. sanguinea") %>%
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
                            group_2=list(species = "F. sanguinea", callow = 1),
                            values_column="mass", group_var="callow")
df$callow <- as.numeric(df$callow)

p1 <- ggplot(aes(x=callow, y=mass, group=callow), data=df)+
  geom_boxplot(col="#1f4a24",
               fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(x=callow+point_x_pos, y=mass, fill=sang_prop), shape=21, cex=3.4, color="black")+
  geom_line(aes(x=callow+point_x_pos, y=mass, group=ID), col="#444444")+
  scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_x_continuous(breaks = c(0,1), labels = c("Callow ants", "Mature ants"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea workers")+
  ylab(expression(paste("Cucitular hydrocarbon mass per individual [",mu,"g]")))+
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
                          values_column="proportion", group_var="callow"
)


df$mature<- as.numeric(df$callow)
p2 <- ggplot(aes(x=as.numeric(callow), y=proportion , group=mature), data=df)+
  geom_boxplot(col="#1f4a24",
               fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(x=as.numeric(callow)+point_x_pos, y=proportion, fill=sang_prop),shape=21,cex=3.4,color="black")+
  #geom_line(aes(x=as.numeric(mature)+pos, y=mass, group=ID), col="#444444")+
  #scale_y_continuous(trans="log2", breaks=2^(0:4))+
  scale_x_continuous(breaks = c(0,1), labels = c("Callow ants", "Mature ants"))+
  scale_fill_gradient("Proportion of\nF. sanguinea ants", low="#555C58", high="#BBEA5E")+
  xlab("F. sanguinea workers")+
  ylab(expression(paste("Cucitular hydrocarbon mass per individual 000 [",mu,"g]")))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=18))+
  theme(legend.title=element_text(size=14),
        legend.text=element_text(size=13))+
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(size = 0.5, linetype = "solid",
                                 colour = "black"))+
  geom_signif(
    y_position = max(df[["proportion"]])*1.1, xmin = 0, xmax = 1,
    annotation = 0.000153, tip_length = 0.02
  )



ggarrange(p1, p2, ncol=2, labels=c("A", "B"))





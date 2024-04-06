library(ggrepel)
library(ggplot2)
library(dplyr)
library(sangchem)
load_globals()
data("development_data")
df <- filter(development_data, caste=="worker",!is.na(head_width), !remarks %in% c("aggression actor", "aggression target", "queenless colony"), species=="F. sanguinea") %>%
  select(colony, census_date, chromatogram_ID,
       callow, mass, head_width, species, sang_prop) %>%
  mutate(mass=mass/(head_width/1300)^2)

df1 <- filter(df, species=="F. sanguinea", callow==1) %>%
  select(colony, census_date, chromatogram_ID, callow, mass, sang_prop) %>%
  group_by(colony,census_date) %>%
  {function(x) {x$pos<-as.factor(group_indices(x)); ungroup(x)}}()

# duplicated samples should converge, thus they should have the same position on x axis
# each level of the factor correspond to the same colony/sampling date combination
levels(df1$pos)<-seq(-0.15,0.15,length.out = length(levels(df1$pos)))
df1$pos <- as.numeric(as.character(df1$pos))
df1$callow_ID <- df1$chromatogram_ID # new column name not to be duplicated by joining with mature ants df


filter(df, species=="F. sanguinea", callow==0) %>%
  select(colony, census_date, chromatogram_ID, callow, mass, sang_prop) -> df2
# pair with callow ants form the same colony and sampled at the same time
df2<-left_join(df2, select(df1, pos, colony, census_date, callow_ID), relationship="many-to-many")
# drop unpaired samples
df2 <- df2[!is.na(df2$callow_ID),]
df2$ID<-1:nrow(df2) # make rows unique before joining to df1

# transer pair ID to callow data frame
df1<-left_join(df1,select(df2, ID,callow_ID))
df<-rbind(df1, df2)
df$mature=df$callow==0

ggplot(aes(x=as.numeric(mature), y=mass , group=mature), data=df)+
  geom_boxplot(col="#1f4a24",
                fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
 geom_point(aes(x=as.numeric(mature)+pos, y=mass, fill=sang_prop),shape=21,cex=3.4,color="black")+
  geom_line(aes(x=as.numeric(mature)+pos, y=mass, group=ID), col="#444444")+
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
                                 colour = "black"))


df$mass <- df$mass*subsample_proportion(peaks_callow, df$chromatogram_ID)


ggplot(aes(x=as.numeric(mature), y=mass , group=mature), data=df2)+
  geom_boxplot(col="#1f4a24",
               fill="#1f4a2466", width=0.6, lwd=1, outlier.shape=NA)+
  geom_point(aes(x=as.numeric(mature)+pos, y=mass, fill=sang_prop),shape=21,cex=3.4,color="black")+
  #geom_line(aes(x=as.numeric(mature)+pos, y=mass, group=ID), col="#444444")+
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
                                 colour = "black"))

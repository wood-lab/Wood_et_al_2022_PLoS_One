#---
title: "Multivariate analysis"
author: "Chelsea Wood"
date: "updated 21 October 2021"
#---


# To assess the response of parasite community composition to host sex and color, 
# we performed a multivariate analysis.


library(vegan)
library(ggplot2)


### Load up your data

rawdata<-read.table("data/lingcod_psite_data.csv",header=T,sep=",")
str(rawdata)


### Extract just the columns that describe parasite abundance

psite_matrix<-rawdata[,-(1:14)]


### Then create a matrix of host factors as well

# One male was identified as "golden brown" G/BR according to Laurel and Bonnie. Reclassify that male as brown.

rawdata$color<-gsub("G/BR","BR",rawdata$color)

n_full <- rawdata %>%
  group_by(sex, color) %>%
  summarise(n = n_distinct(sample_ID))
n_full

rawdata$sex_color_combo<-paste(rawdata$sex,rawdata$color,sep="_")
host_factors<-cbind.data.frame(rawdata$sex_color_combo,rawdata$sex,rawdata$color,rawdata$sampling_location,rawdata$depth_ft)
names(host_factors)<-c("sex_color_combo","sex","color","sampling_location","depth_ft")

adonis(psite_matrix~host_factors$sex_color_combo+host_factors$sampling_location+host_factors$depth_ft,
       permutations=5000, distance="euc")


### Create Supporting Information Figure 1

melted_data$sex_color_combo<-paste(melted_data$sex,melted_data$color,sep="_")

# Tally up averages for each group to make your stacked barplot

stackplot_data <- melted_data %>%
  group_by(sex_color_combo,psite_spp) %>%
  summarise(mean(psite_count))
stackplot_data<-as.data.frame(stackplot_data)
names(stackplot_data)<-c("sex_color_combo","psite_spp","psite_count")

# Plot

ggplot(stackplot_data, aes(fill=psite_spp, y=psite_count, x=sex_color_combo))+
  geom_bar(position="fill", stat="identity", color = "black", size=0.1)+
  scale_color_manual(values=c("black"))+
  ylab("proportion of total parasite community")+
  xlab("sex-color combination")+
  scale_x_discrete(labels=c("F_BL" = "blue females", "F_BR" = "brown females", "M_BL" = "blue males", "M_BR" = "brown males"))+
  theme_minimal()+
  labs(fill="parasite species")+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=14),axis.text.x=element_text(size=14,angle=45,hjust=1),
        axis.title.x=element_text(size=14),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))




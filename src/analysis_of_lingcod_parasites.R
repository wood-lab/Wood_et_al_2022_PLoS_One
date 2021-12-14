#---
title: "Primary analysis of how sex and color influence parasite burden"
author: "Chelsea Wood"
date: "updated 21 October 2021"
#---
  

library(ggplot2)
library(tidyverse)
library(lme4)
library(car)
library(ggeffects)
library(reshape)


### Load up your data

rawdata<-read.table("data/lingcod_psite_data.csv",header=T,sep=",")
str(rawdata)


### Melt this dataset

melted_data<-melt(rawdata,id=c("sample_ID","date_of_collection","sex","TL_cm","fultons_K","HSI","color",
                            "sampling_region","sampling_location","lat","long","depth_ft","date_of_dissection","dissector_initials"))

names(melted_data)<-c("sample_ID","date_of_collection","sex","TL_cm","fultons_K","HSI","color",
                      "sampling_region","sampling_location","lat","long","depth_ft","date_of_dissection","dissector_initials",
                      "psite_spp","psite_count")

str(melted_data)


### How many fish in the full dataset?

n_full <- melted_data %>%
  group_by(sex, color) %>%
  summarise(n = n_distinct(sample_ID))
n_full


### One male was identified as "golden brown" G/BR according to Laurel and Bonnie. Reclassify that male as brown.

melted_data$color<-gsub("G/BR","BR",melted_data$color)

n_full <- melted_data %>%
  group_by(sex, color) %>%
  summarise(n = n_distinct(sample_ID))
n_full


### Re-scale numeric values before analysis

melted_data$psite_count_sc<-scale(melted_data$psite_count)+2
melted_data$lat_sc<-scale(melted_data$lat)+2
melted_data$depth_ft_sc<-scale(melted_data$depth_ft)+2
melted_data$TL_cm_sc<-scale(melted_data$TL_cm)+2


### Use generalized linear mixed models to assess the correlates of parasite burden for 
### the subset of 89 lingcod selected for parasitological analysis.

nbmodel_interaction<-glmer.nb(as.numeric(psite_count_sc)~depth_ft_sc+sex*color+(1|sampling_region/sampling_location/sample_ID)+
                                (1 + sex*color|psite_spp)+offset(log(TL_cm_sc)),data=melted_data,family="nbinom")
summary(nbmodel_interaction)

vif(nbmodel_interaction)


### Create Figure 3b

color_predict<-ggeffect(nbmodel_interaction,c("color","sex"))

color_plot<-ggplot(color_predict,aes(c(1.95,0.95,2.05,1.05),predicted),group=sex,linetype=sex)+
  geom_errorbar(data=color_predict,mapping=aes(x=c(1.95,0.95,2.05,1.05),ymin=conf.low,ymax=conf.high),width=0.03)+
  geom_line(aes(group=group,linetype=group))+
  geom_point(size=6,pch=21,fill=c("cadetblue","burlywood4","cadetblue","burlywood4"))+
  xlab("host coloration")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  labs(linetype="host sex")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color=c("burlywood4","cadetblue")),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=rev(levels(color_predict$x)),labels=c("brown","blue"))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=14))
color_plot


### Create Figure 4

# To do this, extract the random slopes of MPA for each parasite species. 
# This resource was super helpful for figuring this out:
# https://www.r-bloggers.com/getting-the-most-of-mix-models-with-random-slopes/

a<-fixef(nbmodel_interaction)
b<-ranef(nbmodel_interaction,condVar=TRUE)

# Extract the variances of the random effects

qq<-attr(b[[3]],"postVar")
e<-sqrt(qq)
e<-e[4,4,]

# Calculate the CIs

liminf<-(b[[3]][4]+a[5])-(e*2)
mean_<-(b[[3]][4]+a[5])
limsup<-(b[[3]][4]+a[5])+(e*2)

dotchart(mean_[,1],labels=rownames(mean_),cex=0.5)

# add CIs

for(i in 1:nrow(mean_)){
  lines(x=c(liminf[i,1],limsup[i,1]),y=
          c(i,i))
}

raneff_data<-cbind.data.frame(rownames(mean_),mean_,e,(e^2),mean_-(e*2),mean_+(e*2),row.names=NULL)
names(raneff_data)<-c("psite_taxon","interaction","sd","var","min","max")

# Visualize

ranef_plot<-ggplot(raneff_data,aes(psite_taxon,interaction))+
  geom_point(size=3)+
  geom_errorbar(data=raneff_data,mapping=aes(ymin=min,ymax=max),width=0.5)+
  geom_hline(yintercept = -0.71145, linetype = "dotted")+
  coord_flip()+
  xlab("parasite taxon")+
  ylab("interaction term")+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=14),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=rev(c("CHONAR","LEPPRA","LEPBRE","CHALIM","LIRVUL","GNASPP","UDOADU","DERVAR","RHIELO","PROAPE",
                                "LECGIB","BUCEP1","BUCEP2","PODTHE","METACE","TRYSPP","CUCELO","THYADU","NEMLAR","HYSMAG",
                                "NEMWHI","CORWEG","CORCET","ACANT1","ACANT2")),
                   labels=c("CHONAR"=expression(paste(italic("Chondracanthus narium"))),
                            "LEPPRA"=expression(paste(italic("Lepeophtheirus pravipes"))),
                            "LEPBRE"=expression(paste(italic("Lepeophtheirus breviventris"))),
                            "CHALIM"="chalimus copepod larvae",
                            "LIRVUL"=expression(paste(italic("Lironeca vulgaris"))),
                            "GNASPP"="gnathiid spp.",
                            "UDOADU"=expression(paste(italic("Udonella")," spp.")),
                            "DERVAR"=expression(paste(italic("Derogenes varicus"))),
                            "RHIELO"=expression(paste(italic("Rhipidocotyle elongata"))),
                            "PROAPE"=expression(paste(italic("Prosorhynchus apertus"))),
                            "LECGIB"=expression(paste(italic("Lecithaster gibbosus"))),
                            "BUCEP1"="bucephalid sp. 1",
                            "BUCEP2"="bucephalid sp. 2",
                            "PODTHE"=expression(paste(italic("Podocotyle theragrae"))),
                            "METACE"="fin and muscle metacercariae",
                            "TRYSPP"="Trypanorhyncha spp.",
                            "CUCELO"=expression(paste(italic("Cucullanus elongatus"))),
                            "THYADU"=expression(paste(italic("Hysterothylacium aduncum"))),
                            "NEMLAR"="larval nematodes",
                            "HYSMAG"=expression(paste(italic("Hysterothylacium magnum"))),
                            "NEMWHI"="nematode sp. 1",
                            "CORWEG"=expression(paste(italic("Corynosoma wegeneri"))),
                            "CORCET"=expression(paste(italic("Corynosoma cetaceum"))),
                            "ACANT1"="acanthocephalan sp. 1",
                            "ACANT2"="acanthocephalan sp. 2"))
#theme(legend.position="none")
ranef_plot


# Create Supporting Information Table 1, 
# a table summarizing, for each parasite taxon, mean abundance and abundance as a proportion of all parasites.

psite_table <- melted_data %>%
  group_by(psite_spp) %>%
  summarise(sum = sum(psite_count))
View(psite_table)

psite_table_percent <- melted_data %>%
  group_by(psite_spp) %>%
  summarise(percent = 100*(sum(psite_count)/(sum(psite_table$sum))))
View(psite_table_percent)



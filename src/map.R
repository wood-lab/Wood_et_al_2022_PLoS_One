#---
title: "Map"
author: "Chelsea Wood"
date: "updated 21 October 2021"
#---


dev.off()
options(device = "RStudioGD")

library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales) #for transparency
library(rgdal)
library(ggmap)
library(ggsn)
library(cowplot)
library(geosphere)
library(tidyverse)


### Load up your data

rawdata<-read.table("data/lingcod_psite_data.csv",header=T,sep=",")
str(rawdata)


### One male was identified as "golden brown" G/BR according to Laurel and Bonnie. Reclassify that male as brown.

rawdata$color<-gsub("G/BR","BR",rawdata$color)

n_full <- rawdata %>%
  group_by(sex, color) %>%
  summarise(n = n_distinct(sample_ID))


### What is the depth range?

min(rawdata$depth_ft)
max(rawdata$depth_ft)
max(rawdata$depth_ft)-min(rawdata$depth_ft)


### Extract the lats and longs for all of the sites

rawdata$sex_color_combo<-paste(rawdata$sex,rawdata$color,sep="_")
is_M_BL<-vector()

for(i in 1:length(rawdata$sex_color_combo)){
  if(rawdata$sex_color_combo[i] == "M_BL")
    is_M_BL[i]<-"Y"
  else
    is_M_BL[i]<-"N"
}
rawdata$is_M_BL<-is_M_BL

M_BL_depth<-vector()
for(i in 1:length(rawdata$sex_color_combo)){
  if(rawdata$sex_color_combo[i] == "M_BL")
    M_BL_depth[i]<-rawdata$depth_ft[i]
  else
    M_BL_depth[i]<-""
}
rawdata$M_BL_depth<-M_BL_depth


### Replace feet with meters

108/3.28
rawdata$M_BL_depth<-gsub(108,32.9,rawdata$M_BL_depth)
121/3.28
rawdata$M_BL_depth<-gsub(121,36.9,rawdata$M_BL_depth)
35/3.28
rawdata$M_BL_depth<-gsub(35,10.7,rawdata$M_BL_depth)
44/3.28
rawdata$M_BL_depth<-gsub(44,13.4,rawdata$M_BL_depth)

lingcod_sites <- rawdata %>%
  group_by(sampling_region,sampling_location) %>%
  summarise(lat = mean(lat), long = mean(long), n_fish = n_distinct(sample_ID), sex_color_combo = max(is_M_BL), 
            annotation = max(M_BL_depth))


### How many sites?

nlevels(lingcod_sites$sampling_location)


### Now start mapping

bounds<-c(left=-142.5, bottom=30, right=-115 , top=60.5)
full<-get_stamenmap(bounds, zoom=7, maptype = "terrain-background") %>% ggmap()+
  geom_point(data=lingcod_sites, aes(x=lingcod_sites$long,y=lingcod_sites$lat,fill=as.factor(lingcod_sites$sex_color_combo)),size=3,shape=21)+
  scale_fill_manual(values=c("burlywood4","cadetblue1"),guide=FALSE)+
  geom_rect(xmin=-142.5, ymin=53.5, xmax=-130, ymax=60.5, fill="", color="black")+
  geom_rect(xmin=-126, ymin=31.5, xmax=-115, ymax=41.5, fill="", color="black")+
  geom_rect(xmin=-124, ymin=47, xmax=-121.5, ymax=49, fill="", color="black")+
  xlab("longitude (°W)")+
  ylab("latitude (°N)")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16),plot.margin=margin(0,0,0,0))+
  annotate("text", x=c(-127.5,-126.5,-128.5), y=c(57.5,48.25,37), label = c("(b)","(c)","(d)"), size = 8)

bounds<-c(left=-142.5, bottom=53.5, right=-130, top=60.5)
AK_map<-get_stamenmap(bounds, zoom=7, maptype = "terrain-background") %>% ggmap()+
  geom_point(data=lingcod_sites, aes(x=lingcod_sites$long,y=lingcod_sites$lat,fill=as.factor(lingcod_sites$sex_color_combo)),size=3,shape=21,color="black")+
  scale_fill_manual(values=c("burlywood4","cadetblue1"),guide=FALSE)+
  xlab("longitude (°W)")+
  ylab("latitude (°N)")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16))+
  annotate("text", x=(lingcod_sites$long-0.85), y=lingcod_sites$lat, label = lingcod_sites$annotation)

bounds<-c(left=-126, bottom=31.5, right=-115,top=41.5)
CA_map<-get_stamenmap(bounds, zoom=7, maptype = "terrain-background") %>% ggmap()+
  geom_point(data=lingcod_sites, aes(x=lingcod_sites$long,y=lingcod_sites$lat,fill=as.factor(lingcod_sites$sex_color_combo)),size=3,shape=21)+
  scale_fill_manual(values=c("burlywood4","cadetblue1"),guide=FALSE)+
  xlab("longitude (°W)")+
  ylab("latitude (°N)")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16))+
  annotate("text", x=(lingcod_sites$long-0.85), y=lingcod_sites$lat, label = lingcod_sites$annotation)

bounds<-c(left=-124, bottom=47, right=-121.5, top=49)
WA_map<-get_stamenmap(bounds, zoom=7, maptype = "terrain-background") %>% ggmap()+
  geom_point(data=lingcod_sites, aes(x=rev(jitter(lingcod_sites$long,factor=20)),y=rev(jitter(lingcod_sites$lat,factor=20)),
                             fill=rev(as.factor(lingcod_sites$sex_color_combo))),size=3,shape=21,color="black")+
  scale_fill_manual(values=c("burlywood4","cadetblue1"),guide=FALSE)+
  xlab("longitude (°W)")+
  ylab("latitude (°N)")+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=16))+
  annotate("text", x=(lingcod_sites$long-0.21), y=lingcod_sites$lat, label = lingcod_sites$annotation)

final_map<-ggdraw(plot=NULL,xlim=c(0,30),ylim=c(0,21))+
  draw_plot(full,x=0,y=1.5,width=14,height=17)+
  draw_plot(AK_map,x=11,y=10,width=10,height=10)+
  draw_plot(WA_map,x=19,y=10,width=10,height=10)+
  draw_plot(CA_map,x=15,y=0,width=10,height=10)+
  theme(plot.margin=margin(0,0,0,0))+
  draw_text("(a)",size=26,x=3,y=19.5)+
  draw_text("(b)",size=26,x=12.5,y=20)+
  draw_text("(c)",size=26,x=20.5,y=20)+
  draw_text("(d)",size=26,x=16.5,y=9)

final_map



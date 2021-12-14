#---
title: "Analysis of lingcod population (i.e., frequency of sex and color and their influence on body condition)"
author: "Chelsea Wood"
date: "updated 21 October 2021"
#---


library(reshape)
library(lme4)
library(ggeffects)
library(ggplot2)


### Read in data

lingcod_color<-read.table("data/lingcod_population.csv",header=T,sep=",")


### Make a dataset appropriate for chi-squared test

final_data<- lingcod_color %>%
  group_by(Sex,Color) %>%
  summarise(n = n_distinct(Sample.ID))

final_data<-as.data.frame(final_data)

sex<-c("F","M")
blue<-c(265,60)
brown<-c(711,1055)
265+60+711+1055

compiled_data<-cbind.data.frame(sex,blue,brown)
rownames(compiled_data)<-compiled_data$sex
final_data<-compiled_data[,-1]


### Now run your chi-squared test

chisq<-chisq.test(final_data)
chisq


### Create Figure 1a

compiled_data
plot_data<-melt(compiled_data,id="sex")

color_barplot<-ggplot(plot_data,aes(sex,value,fill=variable))+
  geom_bar(position="stack",stat="identity",color="black")+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("host sex")+
  ylab("number of lingcod individuals sampled")+
  theme_minimal()+
  labs(fill="color")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),
        axis.title.y=element_text(size=16),axis.text.x=element_text(size=18),axis.title.x=element_text(size=18),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=14))
color_barplot


### Does Fulton's K differ among blue/brown, male/female lingcod?

str(lingcod_color)
min(lingcod_color$FultonsK,na.rm=T)
max(lingcod_color$FultonsK,na.rm=T)
hist(lingcod_color$FultonsK)

new<-lingcod_color[lingcod_color[,17]<89,]
min(new$FultonsK,na.rm=T)
max(new$FultonsK,na.rm=T)
hist(new$FultonsK)

new$depth_ft_sc<-scale(new$Depth.ft)+2

skinny_model_gaus<-glmer(FultonsK~Sex*Color+depth_ft_sc+(1|Region/Location),data=new,family="gaussian")
summary(skinny_model_gaus)

coefs <- data.frame(coef(summary(skinny_model_gaus)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs


### Create Figure 1c

skinny_predict<-ggeffect(skinny_model_gaus,c("Color","Sex"))

skinny_plot<-ggplot(skinny_predict,aes(c(1.95,0.95,2.05,1.05),predicted),group=sex,linetype=sex)+
  geom_errorbar(data=skinny_predict,mapping=aes(x=c(1.95,0.95,2.05,1.05),ymin=conf.low,ymax=conf.high),width=0.03)+
  geom_line(aes(group=group,linetype=group))+
  geom_point(size=6,pch=21,fill=c("cadetblue","burlywood4","cadetblue","burlywood4"))+
  xlab("host coloration")+
  ylab("Fulton's K")+
  theme_minimal()+
  labs(linetype="host sex")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),
        axis.title.y=element_text(size=18),axis.text.x=element_text(size=18,color=c("burlywood4","cadetblue")),
        axis.title.x=element_text(size=16),panel.background=element_rect(fill="white",color="black"),
        panel.grid.major=element_line(color=NA),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=rev(levels(skinny_predict$x)),labels=c("brown","blue"))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=14))
skinny_plot


### Does hepatosomatic index (HSI) differ among blue/brown, male/female lingcod?

str(lingcod_color)
min(lingcod_color$HSI,na.rm=T)
max(lingcod_color$HSI,na.rm=T)
hist(lingcod_color$HSI)

new<-lingcod_color[lingcod_color[,18]<18,]
new<-new[new[,18]>0,]
min(new$HSI,na.rm=T)
max(new$HSI,na.rm=T)
hist(new$HSI)

new$depth_ft_sc<-scale(new$Depth.ft)+2

new$HSI_trans<-(new$HSI)^(0.25)
hist(new$HSI_trans)

HSI_model_gaus<-glmer(new$HSI_trans~Sex*Color+depth_ft_sc+(1|Region/Location),data=new,family="gaussian")
summary(HSI_model_gaus)

coefs <- data.frame(coef(summary(HSI_model_gaus)))
# use normal distribution to approximate p-value
coefs$p.z <- 2 * (1 - pnorm(abs(coefs$t.value)))
coefs


### Create Figure 1d

HSI_predict<-ggeffect(HSI_model_gaus,c("Color","Sex"))

HSI_plot<-ggplot(HSI_predict,aes(c(1.95,0.95,2.05,1.05),predicted),group=sex,linetype=sex)+
  geom_errorbar(data=HSI_predict,mapping=aes(x=c(1.95,0.95,2.05,1.05),ymin=conf.low,ymax=conf.high),width=0.03)+
  geom_line(aes(group=group,linetype=group))+
  geom_point(size=6,pch=21,fill=c("cadetblue","burlywood4","cadetblue","burlywood4"))+
  xlab("host coloration")+
  ylab("hepatosomatic index")+
  theme_minimal()+
  labs(linetype="host sex")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),
        axis.title.y=element_text(size=18),axis.text.x=element_text(size=18,color=c("burlywood4","cadetblue")),
        axis.title.x=element_text(size=16),panel.background=element_rect(fill="white",color="black"),
        panel.grid.major=element_line(color=NA),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=rev(levels(HSI_predict$x)),labels=c("brown","blue"))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=14))
HSI_plot




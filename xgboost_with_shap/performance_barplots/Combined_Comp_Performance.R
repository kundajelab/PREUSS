rm(list=ls())
data=read.table("Combined_Comp_Performance.tsv",header=TRUE,sep='\t')
library(reshape2)
m=melt(data)
library(ggplot2)
ggplot(data=m,
       aes(x=m$variable,
           y=m$value,
           group=m$Substrate,
           fill=m$Substrate))+
  geom_bar(stat='identity',position='dodge')+
  scale_fill_manual(name="Substrate",values=c("#1b9e77","#d95f02","#7570b3"))+
  xlab("Performance Metric")+
  ylab("Peformance on\nHeld-Out Test Set")+
  theme(axis.text.x = element_text(angle = 90))
  
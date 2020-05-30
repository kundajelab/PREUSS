rm(list=ls())
library(ggplot2)
library(reshape2)
library(preprocessCore)
source("~/helpers.R")
#Quantile normalize TTYH2
ttyh2=read.table("TTYH2-sg3-pairwise-remove-outlier.csv",header=TRUE,sep='\t',row.names=1)
ttyh2_norm=as.data.frame(normalize.quantiles(as.matrix(ttyh2)))
row.names(ttyh2_norm)=row.names(ttyh2)
names(ttyh2_norm)=names(ttyh2)
write.table(ttyh2_norm,file="TTYH2-sg3-pairwise-remove-outlier-qnorm.csv",quote=FALSE,row.names = TRUE,col.names = TRUE)

#Quantile normalize NEIL1
neil1=read.table("NEIL1-ds1-pairwise-remove-outlier.csv",header=TRUE,sep='\t',row.names=1)
neil1_norm=as.data.frame(normalize.quantiles(as.matrix(neil1)))
row.names(neil1_norm)=row.names(neil1)
names(neil1_norm)=names(neil1)
write.table(neil1_norm,file="NEIL1-ds1-pairwise-remove-outlier-qnorm.csv",quote=FALSE,row.names = TRUE,col.names = TRUE)

#check the distribution plots 
m_neil1=melt(neil1)
m_neil1_norm=melt(neil1_norm)
p1=ggplot(data=m_neil1,aes(x=m_neil1$value,
                           group=m_neil1$variable,
                           fill=m_neil1$variable))+
  geom_histogram(position='dodge',bins=100,alpha=0.5)+
  ggtitle("NEIL1")
p2=ggplot(data=m_neil1_norm,aes(x=m_neil1_norm$value,
                                group=m_neil1_norm$variable,
                                fill=m_neil1_norm$variable))+
  geom_histogram(position='dodge',bins=100,alpha=0.5)+
  ggtitle("NEIL1 Quantile Normalized")

m_ttyh2=melt(ttyh2)
m_ttyh2_norm=melt(ttyh2_norm)
p3=ggplot(data=m_ttyh2,aes(x=m_ttyh2$value,
                           group=m_ttyh2$variable,
                           fill=m_ttyh2$variable))+
  geom_histogram(position='dodge',bins=100,alpha=0.5)+
  ggtitle("TTYH2")
p4=ggplot(data=m_ttyh2_norm,aes(x=m_ttyh2_norm$value,
                                group=m_ttyh2_norm$variable,
                                fill=m_ttyh2_norm$variable))+
  geom_histogram(position='dodge',bins=100,alpha=0.5)+
  ggtitle("TTYH2 Quantile Normalized")

multiplot(p1,p2,p3,p4,cols=2)
rm(list=ls())
library(ggplot2)
data=read.table("FeatureGroupContributionsForPlot.csv",header=TRUE,sep='\t')
data$Subset=factor(data$Subset,level=c("Motif","Downstream", "Upstream","EditingSiteStructure","Thermodynamics","StructureSequence"))
ggplot(data=data,aes(x=data$Subset,y=data$Contribution,size=data$Contribution,color=data$Substrate))+
  geom_point()+
  scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_color_manual(values=c('#1b9e77','#d95f02','#7570b3'))+
  scale_size(range = c(1, 12))+
  coord_flip()+
  theme_bw()
  

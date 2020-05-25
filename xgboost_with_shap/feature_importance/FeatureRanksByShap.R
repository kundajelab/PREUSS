rm(list=ls())
library(ggplot2)
data=read.table("FeatureRanksByShap.tsv",header=TRUE,sep='\t')
data=data[order(data$Rank),]
#take the top 50 
data=data[data$Rank<21,]
data$Feature=factor(data$Feature,levels=data$Feature)
pdf("FeatureRanksByMeanAbsSHAP.pdf",height=6,width=15)
print(ggplot(data=data,aes(x=data$Feature,y=data$MeanImpact))+
  geom_bar(stat='identity',color='black',position=position_dodge())+
  geom_errorbar(aes(ymin=data$MeanImpact-data$StdImpact,
                    ymax=data$MeanImpact+data$StdImpact),
                width=0.2,
                position=position_dodge(0.9))+
  geom_point(aes(x=data$Feature,y=data$AJUBA,color="AJUBA"))+
  geom_point(aes(x=data$Feature,y=data$NEIL1,color="NEIL1"))+
  geom_point(aes(x=data$Feature,y=data$TTYH2,color="TTYH2"))+
  scale_color_manual(name="Substrate", values=c("#1b9e77","#d95f02","#7570b3"))+
  labs(title="Feature Ranks by mean(|SHAP|)",
       x="Feature Label",
       y = "% Contribution to SHAP values")+
  theme_classic(20)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off() 


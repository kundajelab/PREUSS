rm(list=ls())
library(ggplot2)
data=read.table("FeatureRanksByShap.tsv",header=TRUE,sep='\t')
#take the top 50 
data=data[data$MeanRank<51,]
data$Feature=factor(data$Feature,levels=data$Feature)
ggplot(data=data,aes(x=data$Feature,y=data$Mean_.ContribToPrediction))+
  geom_bar(stat='identity',color='black',position=position_dodge())+
  geom_errorbar(aes(ymin=data$Mean_.ContribToPrediction-data$Std_.ContribToPrediction,
                    ymax=data$Mean_.ContribToPrediction+data$Std_.ContribToPrediction),
                width=0.2,
                position=position_dodge(0.9))+
  geom_point(aes(x=data$Feature,y=data$AJUBA_.ContribToPrediction,color="AJUBA"))+
  geom_point(aes(x=data$Feature,y=data$NEIL1_.ContribToPrediction,color="NEIL1"))+
  geom_point(aes(x=data$Feature,y=data$TTYH2_.ContribToPrediction,color="TTYH2"))+
  scale_color_manual(name="Substrate", values=c("#1b9e77","#d95f02","#7570b3"))+
  labs(title="Feature Ranks by mean(|SHAP|)",
       x="Feature Label",
       y = "% Contribution to SHAP values")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylim(0,50)


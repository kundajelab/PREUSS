rm(list=ls())
library(ggplot2)
data=read.table("feature_ranks.tsv",header=TRUE,sep='\t')

#sort by mean rank 
mean_rank=data[order(data$Mean_Rank),]
mean_rank$Feature=factor(mean_rank$Feature,levels=mean_rank$Feature)

p1=ggplot(data=mean_rank,
       aes(x=mean_rank$Feature,
           y=mean_rank$Mean_Rank))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean_rank$Mean_Rank-sqrt(mean_rank$Var_Rank), ymax=mean_rank$Mean_Rank+sqrt(mean_rank$Var_Rank)), width=.2,
                position=position_dodge(.5))+
  geom_point(aes(x=mean_rank$Feature,
               y=mean_rank$NEIL1_Rank,
               color="NEIL1"))+
  geom_point(aes(x=mean_rank$Feature,
               y=mean_rank$AJUBA_Rank,
               color="AJUBA"))+
  geom_point(aes(x=mean_rank$Feature,
               y=mean_rank$TTYH2_Rank,
               color="TTYH2"))+
  xlab("Feature")+
  ylab("Mean Rank +/-\nStd. Dev.")+
  theme_bw(20)+
  theme(axis.text.x = element_text(angle=90, hjust=1))

#sort by mean increase in MSE 
mean_mse=data[order(data$Mean_.IncMSE,decreasing = TRUE),]
mean_mse$Feature=factor(mean_mse$Feature,levels=mean_mse$Feature)

p2=ggplot(data=mean_mse,
          aes(x=mean_mse$Feature,
              y=mean_mse$Mean_.IncMSE))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=mean_mse$Mean_.IncMSE-sqrt(mean_mse$Var_.IncMSE),
                    ymax=mean_mse$Mean_.IncMSE+sqrt(mean_mse$Var_.IncMSE)), 
                width=.2,
                position=position_dodge(.5))+
  geom_point(aes(x=mean_mse$Feature,
                 y=mean_mse$NEIL1_.IncMSE,
                 color="NEIL1"))+
  geom_point(aes(x=mean_mse$Feature,
                 y=mean_mse$AJUBA_.IncMSE,
                 color="AJUBA"))+
  geom_point(aes(x=mean_mse$Feature,
                 y=mean_mse$TTYH2_.IncMSE,
                 color="TTYH2"))+
  xlab("Feature")+
  ylab("Mean %Inc MSE +/-\nStd. Dev.")+
  theme_bw(20)+
  theme(axis.text.x = element_text(angle=90, hjust=1))

svg("FeaturesSortedByRank.svg",height=8,width=15)
print(p1)
dev.off() 

svg("FeaturesSortedByMSE.svg",height=8,width=15)
print(p2)
dev.off() 

#sort by rank variance 
var_rank=data[order(data$Var_Rank),]
var_rank$Feature=factor(var_rank$Feature,levels=var_rank$Feature)

p3=ggplot(data=var_rank,
          aes(x=var_rank$Feature,
              y=sqrt(var_rank$Var_Rank)))+
  geom_bar(stat="identity")+
  geom_point(aes(x=var_rank$Feature,
                 y=var_rank$NEIL1_Rank,
                 color="NEIL1"))+
  geom_point(aes(x=var_rank$Feature,
                 y=var_rank$AJUBA_Rank,
                 color="AJUBA"))+
  geom_point(aes(x=var_rank$Feature,
                 y=var_rank$TTYH2_Rank,
                 color="TTYH2"))+
  xlab("Feature")+
  ylab("Rank Standard Error")+
  theme_bw(20)+
  theme(axis.text.x = element_text(angle=90, hjust=1))

#sort by MSE variance 
var_mse=data[order(data$Var_.IncMSE),]
var_mse$Feature=factor(var_mse$Feature,levels=var_mse$Feature)

p4=ggplot(data=var_mse,
          aes(x=var_mse$Feature,
              y=sqrt(var_mse$Var_.IncMSE)))+
  geom_bar(stat="identity")+
  geom_point(aes(x=var_mse$Feature,
                 y=var_mse$NEIL1_.IncMSE,
                 color="NEIL1"))+
  geom_point(aes(x=var_mse$Feature,
                 y=var_mse$AJUBA_.IncMSE,
                 color="AJUBA"))+
  geom_point(aes(x=var_mse$Feature,
                 y=var_mse$TTYH2_.IncMSE,
                 color="TTYH2"))+
  xlab("Feature")+
  ylab("%IncMSE Standard Error")+
  theme_bw(20)+
  theme(axis.text.x = element_text(angle=90, hjust=1))


svg("FeaturesSortedByVarianceOfRank.svg",height=8,width=15)
print(p3)
dev.off() 

svg("FeaturesSortedByVarianceOfMSE.svg",height=8,width=15)
print(p4)
dev.off() 
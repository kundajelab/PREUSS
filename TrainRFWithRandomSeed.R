#!/usr/bin/env Rscript

library(randomForest)
library(gbm)
library(reshape2)
library(ggplot2)
library(optparse)
options(warn=-1)
##define function to calculate mean absolute percent error 
mape <- function(actual,pred){
  #pseudocount 
  actual=actual+0.001
  pred=pred+0.001
  mape <- mean(abs((actual - pred)/actual))*100
  print(paste("MAPE:",mape))
  return (mape)
}

medape <- function(actual,pred){
  #pseudocount 
  actual=actual+0.001
  pred=pred+0.001
  medape=median(abs((actual-pred)/actual))*100
  print(paste("MED APE: ",medape))
  return (medape)
}

## Read in user arguments ## 
option_list=list(
  make_option(c("--features_mat",type="character",default=NULL,help="Path to feature matrix")),
  make_option(c("--output_dir",type="character",default=".",help="Directory to store output files")),
  make_option(c("--n_iter",type="integer",default=10,help="Number of train/test splits to generate for training")),
  make_option(c("--seed",type="integer",default=1234,help="Random seed integer")),
  make_option(c("--train_split_size",type="double",default=as.double(0.8),help="Fraction of feature matrix entries to use for training")),
  make_option(c("--number_trees",type="integer",default=500,help="Number of trees to grow in random forest")),
  make_option(c("--n_features_to_plot",type="integer",default=10,help="N most important features to plot"))
  
  
)
opt_parser=OptionParser(option_list=option_list)
opt=parse_args(opt_parser)


if (is.null(opt$features_mat)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
if(is.null(opt$output_dir)){
  output_dir='.'
}else{
  output_dir=opt$output_dir
}

if(is.null(opt$n_iter)){
  n_iter=10
}else{
  n_iter=as.integer(opt$n_iter)
}
if(is.null(opt$seed)){
  seed_int=1234
}else{
  seed_int=as.integer(opt$seed)
}
if (is.null(opt$train_split_size)){
  train_split_size=0.8
}else{
  train_split_size=as.double(opt$train_split_size)
}
test_split_size=1-train_split_size

if (is.null(opt$number_trees)){
  number_trees=500
}else{
  number_trees=as.integer(opt$number_trees)
}
if (is.null(opt$n_features_to_plot)){
  n_features_to_plot=10
}else{
  n_features_to_plot=as.integer(opt$n_features_to_plot)
}
features_mat=opt$features_mat



##create the output directory if it doesn't exist 
dir.create(paste(getwd(),output_dir,sep='/'),showWarnings = FALSE)

##Load Feature Matrix 
data=read.table(features_mat,header=TRUE,sep='\t',row.names = 1)
data=data[is.na(data$editing_value)==FALSE,]

#set random seed 
set.seed(seed_int)
print(paste("seed:",seed_int))

#drop any columns that are composed of NA entirely
data=data[colSums(is.na(data)) < nrow(data)]

#fill in missing values via imputation 
data <- na.roughfix(data)

#Remove any columns with 0 variance 
data=data[sapply(data,var)>0]


## Train and Test Splits

#number of iterations of data splitting + random forest 
#define train/test split sizes and number of iterations. 
print(paste("n_iter:",n_iter))
print(paste("train_split_size:",train_split_size))
print(paste("test_split_size:",test_split_size))

n_train=floor(train_split_size*nrow(data))
n_test=nrow(data)-n_train
## Keep track of training outputs 
train_truth=matrix(0,nrow=n_train,ncol=n_iter)
test_truth=matrix(0,nrow=n_test,ncol=n_iter) 

train_pred=matrix(0,nrow=n_train,ncol=n_iter)
test_pred=matrix(0,nrow=n_test,ncol=n_iter) 

train_mse=matrix(0,nrow=n_iter,ncol=1)
test_mse=matrix(0,nrow=n_iter,ncol=1)

train_mape=matrix(0,nrow=n_iter,ncol=1)
test_mape=matrix(0,nrow=n_iter,ncol=1)

train_medape=matrix(0,nrow=n_iter,ncol=1)
test_medape=matrix(0,nrow=n_iter,ncol=1)

train_var_explained=matrix(0,nrow=n_iter,ncol=1) 
test_var_explained=matrix(0,nrow=n_iter,ncol=1) 

forests=list()
## Train the ensemble model 
for(iter in seq(1,n_iter)){
  
  train_indices=sample(nrow(data),n_train,replace=FALSE)  
  train_split=data[train_indices,]
  test_split=data[-train_indices,]
  
  ytrain=train_split$editing_value
  xtrain=train_split[,2:ncol(train_split)]
  ytest=test_split$editing_value
  xtest=test_split[,2:ncol(test_split)]
  
  forest=randomForest(y=ytrain,
                      x=xtrain,
                      xtest=xtest,
                      ytest=ytest,
                      keep.forest=TRUE,
                      importance=TRUE,
                      ntree=number_trees)
  
  #extract predictions and performance metrics from the random forest 
  predictions_training_data=forest$predicted
  predictions_test_data=forest$test$predicted
  feat_importance=importance(forest)
  train_mape_iter=mape(ytrain,predictions_training_data)
  train_medape_iter=medape(ytrain,predictions_training_data)
  train_mse_iter=mean(forest$mse)
  
  train_var_explained_iter=mean(forest$rsq)
  test_mape_iter=mape(ytest,predictions_test_data)
  test_medape_iter=medape(ytest,predictions_test_data)
  test_mse_iter=mean(forest$test$mse)
  test_var_explained_iter=mean(forest$test$rsq)
  
  #append all metrics to running list 
  train_truth[,iter]=ytrain
  test_truth[,iter]=ytest 
  train_pred[,iter]=predictions_training_data
  test_pred[,iter]=predictions_test_data
  train_mape[iter,]=train_mape_iter
  train_medape[iter,]=train_medape_iter
  train_mse[iter,]=train_mse_iter 
  test_mape[iter,]=test_mape_iter
  test_medape[iter,]=test_medape_iter
  test_mse[iter,]=test_mse_iter
  train_var_explained[iter,]=train_var_explained_iter
  test_var_explained[iter,]=test_var_explained_iter
  forests[[as.character(iter)]]=forest
}

## Quantify the MSE 
print(paste("Training MSE (Ensembled):",mean(train_mse)))
print(paste("Test MSE (Ensembled):",mean(test_mse)))
print(paste("Training MAPE (Ensembled):",mean(train_mape)))
print(paste("Test MAPE (Ensembled):",mean(test_mape)))
print(paste("Training MEDIAN APE (Ensembled):",mean(train_medape)))
print(paste("Test MEDIAN APE (Ensembled):",mean(test_medape)))

mse_df=data.frame(train_mse,test_mse)
mape_df=data.frame(train_mape,test_mape)
medape_df=data.frame(train_medape,test_medape)

names(mse_df)=c("Training MSE","Test MSE")
names(mape_df)=c("Training MAPE","Test MAPE")
names(medape_df)=c("Training MED.APE","Test MED.APE")

mse_df=melt(mse_df)
p1=ggplot(data=mse_df,
          aes(x=mse_df$variable,y=mse_df$value))+
  geom_boxplot()+
  xlab("Data Split")+
  ylab("MSE (Ensembled)")+
  ggtitle(paste("Training MSE:",round(mean(train_mse),3),"; Test MSE:",round(mean(test_mse),3)))
print(p1)
svg(paste(output_dir,"MSE.svg",sep='/'),height=6,width=6)
print(p1)
dev.off()

mape_df=melt(mape_df)
p2=ggplot(data=mape_df,
          aes(x=mape_df$variable,y=mape_df$value))+
  geom_boxplot()+
  xlab("Data Split")+
  ylab("MAPE (Ensembled)")+
  ggtitle(paste("Training MAPE:",round(mean(train_mape),3),"; Test MAPE:",round(mean(test_mape),3)))
print(p2)
svg(paste(output_dir,"MAPE.svg",sep='/'),height=6,width=6)
print(p2)
dev.off()

medape_df=melt(medape_df)
p21=ggplot(data=medape_df,
          aes(x=medape_df$variable,y=medape_df$value))+
  geom_boxplot()+
  xlab("Data Split")+
  ylab("MEDIAN APE (Ensembled)")+
  ggtitle(paste("Training MEDIAN APE:",round(mean(train_medape),3),"; Test MAPE:",round(mean(test_medape),3)))
print(p21)
svg(paste(output_dir,"MEDAPE.svg",sep='/'),height=6,width=6)
print(p21)
dev.off()

## Quantify % of Variance Explained in the Data by the Random Forest 
print(paste("Training Var Explained (Ensembled):",mean(train_var_explained)))
print(paste("Test Var Explained (Ensembled):",mean(test_var_explained)))
var_df=data.frame(train_var_explained,test_var_explained)
names(var_df)=c("Training Var Explained","Test Var Explained")
var_df=melt(var_df)
p3=ggplot(data=var_df,
          aes(x=var_df$variable,y=var_df$value))+
  geom_boxplot()+
  xlab("Data Split")+
  ylab("% Variance Explained (Ensembled)")+
  ggtitle(paste("Training Var Explained:",round(mean(train_var_explained),3),"; Test Var Explained:",round(mean(test_var_explained),3)))
print(p3)
svg(paste(output_dir,"VarianceExplained.svg",sep='/'),height=6,width=6)
print(p3)
dev.off()

## Train and Test Split Predictions 
#color-code by source of data 
source("helpers.R")
v_train_truth=as.vector(train_truth)
v_train_pred=as.vector(train_pred)
v_test_truth=as.vector(test_truth)
v_test_pred=as.vector(test_pred)

v_train=data.frame(v_train_truth,v_train_pred)
names(v_train)=c("Observed","Predicted")
v_test=data.frame(v_test_truth,v_test_pred)
names(v_test)=c("Observed","Predicted")

p4=ggplot(data=v_train,
          aes(x=v_train$Observed,
              y=v_train$Predicted))+
  geom_point(alpha=0.3)+
  geom_abline()+
  xlab("Observed")+
  ylab("Predicted")+
  ggtitle(paste("Training Splits Editing Levels (Ensembled);",
                "\nMSE:",
                round(mean(train_mse),3),
                "\nMAPE:",
                round(mean(train_mape),3),
                "\nMED APE:",
                round(mean(train_medape),3),
                "; \nVar Explained:",
                round(mean(train_var_explained),3)))+
  theme_bw(10)

p5=ggplot(data=v_test,
          aes(x=v_test$Observed,
              y=v_test$Predicted))+
  geom_point(alpha=0.3)+
  geom_abline()+
  xlab("Observed")+
  ylab("Predicted")+
  ggtitle(paste("Test Splits Editing Levels (Ensembled);",
                "\nMSE:",
                round(mean(test_mse),3),
                "\nMAPE:",
                round(mean(test_mape),3),
                "\nMEDAPE:",
                round(mean(test_medape),3),
                "; \nVar Explained:",
                round(mean(test_var_explained),3)))+
  theme_bw(10)
multiplot(p4,p5,cols=2)

svg(paste(output_dir,"TrainingSplitPerformance.svg",sep='/'),width=6,height=6)
print(p4)
dev.off() 

svg(paste(output_dir,"TestSplitPerformance.svg",sep='/'),width=6,height=6)
print(p5)
dev.off() 


## Write Truth and Predicted Labels to TSV 
write.table(v_train,file=paste(output_dir,"TrainingDataPerformance.tsv",sep='/'),col.names = TRUE,row.names=FALSE,sep='\t',quote=FALSE)
write.table(v_test,file=paste(output_dir,"TestDataPerformance.tsv",sep='/'),col.names = TRUE,row.names=FALSE,sep='\t',quote=FALSE)

## Merge the forests into one ensembl to calculate feature importance 
all_forests=do.call(combine,forests)
combined_importance=importance(all_forests)
print(combined_importance)
#varImpPlot(all_forests,n.var=10,main="Feature Rank (top 10) from Ensemble of Random Forests")
write.table(combined_importance,file=paste(output_dir,"FeatureImportance.tsv",sep='/'),sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)

## Generate Feature vs Editing Level plots for Top 10 Features 
for(feature in row.names(importance)[1:n_features_to_plot])
{
  feat_vals=data[feature]
  editing_vals=data$editing_value
  feat_df=data.frame(feat_vals,editing_vals)
  names(feat_df)=c("FeatureValue","EditingLevel")
  p6=ggplot(data=feat_df,
            aes(x=feat_df$FeatureValue,
                y=feat_df$EditingLevel))+
    geom_point()+
    xlab(as.character(feature))+
    ylab("Editing Level")
  print(p6)
  outputfname=paste(output_dir,feature,sep='/')
  outputfname=paste(outputfname,'svg',sep='.')
  svg(outputfname,height=6,width=6)
  print(p6)
  dev.off() 
}

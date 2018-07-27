rm(list=ls())
library(ggplot2)
library(randomForest)
#load the feature matrix 
data=read.table("rf_features.txt",header=TRUE,sep='\t')
#remove values with no NaN in editing level
data=data[is.nan(data$ave_editing_level)==FALSE,]

#get training and test splits 
set.seed(1)
train_indices=sample(nrow(data),0.8*nrow(data),replace=FALSE)
train_split=data[train_indices,]
test_split=data[-train_indices,]

#train rf 
bag.data=randomForest(ave_editing_level~.,data=data,subset=train_indices,mtry=19,importance=TRUE)
print(bag.data)

#get predictions on training & test data 
yhat.bag.train=predict(bag.data,newdata=train_split)
yhat.bag.test=predict(bag.data,newdata=test_split)
plot(yhat.bag.train,train_split$ave_editing_level)
abline(0,1)
mean((yhat.bag.train-train_split$ave_editing_level)^2)

plot(yhat.bag.test,test_split$ave_editing_level)
abline(0,1)
mean((yhat.bag.test-test_split$ave_editing_level)^2)

#get the feature importance 
importance (rf.data )
varImpPlot(rf.data)
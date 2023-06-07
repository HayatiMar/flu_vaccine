library(e1071)
library(DMwR)
library(outliers)
library(ggplot2)
library("ROCR")
library(caret)
library("gbm")
library("Boruta")
source("code/Tree_Statistics.R")
source("code/general_functions.R")
#data containing tips before 28/2/2017

data = read.csv("2020/mycurrentdata2020_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
#length(which(data$outcome==1))/nrow(data)
Aux_data = read.csv("2020/Aux_dataNA2020.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#we remove the very past data since we are not sure about their labels!
verypast_ind = which(as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) <="2010-1-1")
verypast_name = Aux_data[verypast_ind ,1]
ii = match(verypast_name,data$tip)
length(which(!is.na(ii)))
verypast_ii = ii[which(!is.na(ii))]
#data$outcome[verypast_ii]
length(which(data$outcome[verypast_ii]==1))
length(verypast_ii)
data = data[-verypast_ii,]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#keep the unique rows
data_tip = data$tip
data = data[,-c(80)]
dup_ind = which(duplicated(data))
length(dup_ind)
data = unique(data)
data_tip = data_tip[-dup_ind]
length(data_tip)
length(which(data$outcome==1))/length(data$outcome)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#remove the outliers
set.seed(17)

alloutliers = c(842,1133,880,1072, 306,14 ,115 ,104 ,108,118 , 102 , 110 ,111,107,13015 ,339,341,485,486, 15022, 15033,
                18204, 18584,5133,4128, 21863, 21124,  8611,  6178 ,18498, 18455, 18595, 18577 ,18583 ,18800 ,21125, 21122 ,21121 ,26015, 21507 ,22478, 17370,
                17369 ,26396, 13019 ,13020, 13038, 18684, 18742 ,21758)

outliers = alloutliers[1:25]
print(outliers)
data=data[-outliers,]
data_tip = data_tip[-outliers]
length(data_tip)
which(is.na(data),arr.ind = TRUE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#divide the data to test and train -> test on a clade (76876) and train on other clades
load("2020/flutree2020-2.Rdata")
test_subtree = extract.clade(tree,76876)
test_subtree
length(test_subtree$tip.label)
test_names = test_subtree$tip.label
test_ind = match(test_names,data_tip)
#some of the tips are after Feb 2017, so they are in the prediction table. And some of them are those that we do not have their
#NA sequesne
length(which(is.na(test_ind)))
testi = test_ind[which(!is.na(test_ind))]
length(testi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ind = which(data$outcome==1)
data$outcome[ind] = 0
data$outcome[-ind] =1
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
length(which(data$outcome==0))/length(data$outcome)
data = scaleData(data)
test = data[testi,]
test = test[sample(seq_len(nrow(test)), size = nrow(test)),]
train = data[-testi,]
train = train[sample(seq_len(nrow(train)), size = nrow(train)),]
which(is.na(train),arr.ind = TRUE)
#train = train[-unique(which(is.na(train),arr.ind = TRUE)[,1]),]
which(is.na(test),arr.ind = TRUE)
length(which(train$outcome==0))/length(train$outcome)
length(which(test$outcome==0))/length(test$outcome)
nrow(test)/(nrow(train)+nrow(test))
#==============================================================================================
#experiments with all the features
svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.005,cost =128,degree=5,coef0=0)
print(svm.fit)
mypred=predict(svm.fit,test)
table(mypred, test$outcome) 
sum(mypred==test$outcome)/nrow(test)
confusionMatrix(mypred, test$outcome) 
svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
svmmodel.probs<-attr(svmmodel.predict,"decision.values")
svmmodel.class<-predict(svm.fit,test,type="class")
svmmodel.outcome<-test$outcome
#roc analysis for test data
svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.outcome)
svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
round(svmmodel.auc,2)
plot(svmmodel.performance,type="l", lwd=1,col="red",cex.lab=1, yaxt="n", xaxt="n", ylab="",xlab="")
#title("ROC curve for classificaton using SVM \n with different kernels" , adj = 0.5, line = 1,cex.main=1)
#axis(2,cex.axis=0.9)
mtext("True Positive rate", side=2, line=2.3, cex=1)
mtext("False Positive rate", side=1, line=2.5, cex=1)
abline(0,1,col="grey82",  lwd=0.8)
par(new=TRUE)
#==============================================================================================
#==============================================================================================
#==============================================================================================
#experiments with all the features except NA
data = read.csv("2020/mycurrentdata2020_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
data |>
  select(c(-ends_with("_NA"),"outcome")) -> data

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#we remove the very past data since we are not sure about their labels!
verypast_ind = which(as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) <="2010-1-1")
verypast_name = Aux_data[verypast_ind ,1]
ii = match(verypast_name,data$tip)
length(which(!is.na(ii)))
verypast_ii = ii[which(!is.na(ii))]
#data$outcome[verypast_ii]
length(which(data$outcome[verypast_ii]==1))
length(verypast_ii)
data = data[-verypast_ii,]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#keep the unique rows
data_tip = data$tip
data = data[,-c(41)]
dup_ind = which(duplicated(data))
length(dup_ind)
data = unique(data)
data_tip = data_tip[-dup_ind]
length(data_tip)
length(which(data$outcome==1))/length(data$outcome)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#remove the outliers
set.seed(17)

alloutliers = c(842,1133,880,1072, 306,14 ,115 ,104 ,108,118 , 102 , 110 ,111,107,13015 ,339,341,485,486, 15022, 15033,
                18204, 18584,5133,4128, 21863, 21124,  8611,  6178 ,18498, 18455, 18595, 18577 ,18583 ,18800 ,21125, 21122 ,21121 ,26015, 21507 ,22478, 17370,
                17369 ,26396, 13019 ,13020, 13038, 18684, 18742 ,21758)

outliers = alloutliers[1:25]
print(outliers)
data=data[-outliers,]
data_tip = data_tip[-outliers]
length(data_tip)
which(is.na(data),arr.ind = TRUE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#divide the data to test and train -> test on a clade (76876) and train on other clades
load("2020/flutree2020-2.Rdata")
test_subtree = extract.clade(tree,76876)
test_subtree
length(test_subtree$tip.label)
test_names = test_subtree$tip.label
test_ind = match(test_names,data_tip)
#some of the tips are after Feb 2017, so they are in the prediction table. And some of them are those that we do not have their
#NA sequesne
length(which(is.na(test_ind)))
testi = test_ind[which(!is.na(test_ind))]
length(testi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ind = which(data$outcome==1)
data$outcome[ind] = 0
data$outcome[-ind] =1
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
length(which(data$outcome==0))/length(data$outcome)
data = scaleData(data)
test = data[testi,]
test = test[sample(seq_len(nrow(test)), size = nrow(test)),]
train = data[-testi,]
train = train[sample(seq_len(nrow(train)), size = nrow(train)),]
which(is.na(train),arr.ind = TRUE)
#train = train[-unique(which(is.na(train),arr.ind = TRUE)[,1]),]
which(is.na(test),arr.ind = TRUE)
length(which(train$outcome==0))/length(train$outcome)
length(which(test$outcome==0))/length(test$outcome)
nrow(test)/(nrow(train)+nrow(test))

svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.005,cost =128,degree=5,coef0=0)
print(svm.fit)
mypred=predict(svm.fit,test)
table(mypred, test$outcome) 
sum(mypred==test$outcome)/nrow(test)
confusionMatrix(mypred, test$outcome) 
svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
svmmodel.probs<-attr(svmmodel.predict,"decision.values")
svmmodel.class<-predict(svm.fit,test,type="class")
svmmodel.outcome<-test$outcome
#roc analysis for test data
svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.outcome)
svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
round(svmmodel.auc,2)
plot(svmmodel.performance,type="l", lwd=1,col="blue",cex.lab=1, yaxt="n", xaxt="n", ylab="",xlab="")
#title("ROC curve for classificaton using SVM \n with different kernels" , adj = 0.5, line = 1,cex.main=1)
#axis(2,cex.axis=0.9)
mtext("True Positive rate", side=2, line=2.3, cex=1)
mtext("False Positive rate", side=1, line=2.5, cex=1)
abline(0,1,col="grey82",  lwd=0.8)
par(new=TRUE)

#==============================================================================================
#==============================================================================================
#==============================================================================================
#experiments with all the features except epitope
data = read.csv("2020/mycurrentdata2020_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
data |>
  select(c(-matches("epitope"),"outcome")) -> data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#we remove the very past data since we are not sure about their labels!
verypast_ind = which(as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) <="2010-1-1")
verypast_name = Aux_data[verypast_ind ,1]
ii = match(verypast_name,data$tip)
length(which(!is.na(ii)))
verypast_ii = ii[which(!is.na(ii))]
#data$outcome[verypast_ii]
length(which(data$outcome[verypast_ii]==1))
length(verypast_ii)
data = data[-verypast_ii,]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#keep the unique rows
data_tip = data$tip
data = data[,-79]
dup_ind = which(duplicated(data))
length(dup_ind)
data = unique(data)
data_tip = data_tip[-dup_ind]
length(data_tip)
length(which(data$outcome==1))/length(data$outcome)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#remove the outliers
set.seed(17)

alloutliers = c(842,1133,880,1072, 306,14 ,115 ,104 ,108,118 , 102 , 110 ,111,107,13015 ,339,341,485,486, 15022, 15033,
                18204, 18584,5133,4128, 21863, 21124,  8611,  6178 ,18498, 18455, 18595, 18577 ,18583 ,18800 ,21125, 21122 ,21121 ,26015, 21507 ,22478, 17370,
                17369 ,26396, 13019 ,13020, 13038, 18684, 18742 ,21758)

outliers = alloutliers[1:25]
print(outliers)
data=data[-outliers,]
data_tip = data_tip[-outliers]
length(data_tip)
which(is.na(data),arr.ind = TRUE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#divide the data to test and train -> test on a clade (76876) and train on other clades
load("2020/flutree2020-2.Rdata")
test_subtree = extract.clade(tree,76876)
test_subtree
length(test_subtree$tip.label)
test_names = test_subtree$tip.label
test_ind = match(test_names,data_tip)
#some of the tips are after Feb 2017, so they are in the prediction table. And some of them are those that we do not have their
#NA sequesne
length(which(is.na(test_ind)))
testi = test_ind[which(!is.na(test_ind))]
length(testi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ind = which(data$outcome==1)
data$outcome[ind] = 0
data$outcome[-ind] =1
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
length(which(data$outcome==0))/length(data$outcome)
data = scaleData(data)
test = data[testi,]
test = test[sample(seq_len(nrow(test)), size = nrow(test)),]
train = data[-testi,]
train = train[sample(seq_len(nrow(train)), size = nrow(train)),]
which(is.na(train),arr.ind = TRUE)
#train = train[-unique(which(is.na(train),arr.ind = TRUE)[,1]),]
which(is.na(test),arr.ind = TRUE)
length(which(train$outcome==0))/length(train$outcome)
length(which(test$outcome==0))/length(test$outcome)
nrow(test)/(nrow(train)+nrow(test))

svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.005,cost =128,degree=5,coef0=0)
print(svm.fit)
mypred=predict(svm.fit,test)
table(mypred, test$outcome) 
sum(mypred==test$outcome)/nrow(test)
confusionMatrix(mypred, test$outcome) 
svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
svmmodel.probs<-attr(svmmodel.predict,"decision.values")
svmmodel.class<-predict(svm.fit,test,type="class")
svmmodel.outcome<-test$outcome
#roc analysis for test data
svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.outcome)
svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
round(svmmodel.auc,2)

plot(svmmodel.performance,type="l", lwd=1,col="green",cex.lab=1, yaxt="n", xaxt="n", ylab="",xlab="")
title("ROC curve for classification using \n different sets of features " , adj = 0.5, line = 1,cex.main=0.8)
#axis(2,cex.axis=0.9)
mtext("True Positive rate", side=2, line=2.3, cex=1)
mtext("False Positive rate", side=1, line=2.5, cex=1)
abline(0,1,col="grey82",  lwd=0.8)
par(new=TRUE)
legend("bottomleft", legend=c("All features: AUC = 0.89-Accuracy:89%","All features without NA: AUC = 0.72-Accuracy:77%","All feature without epitope: AUC = 0.92-Accuracy:89%"), fill=c("red","blue","green"), cex=0.7)


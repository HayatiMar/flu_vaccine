library(e1071)
library(ape)
library(phangorn)
library(apTreeshape)
library(DMwR2)
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
length(which(data$outcome==1))/nrow(data)
Aux_data = read.csv("2020/Aux_dataNA2020.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
#data[outliers,]
data=data[-outliers,]
data_tip = data_tip[-outliers]
length(data_tip)
which(is.na(data),arr.ind = TRUE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#divide the data to test -> test on a clade (102376,76876) and train on other clades
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
#=========================================================================================
#Random Forest
set.seed(17)
library(randomForest)
# Perform training:
rf.fit = randomForest(outcome ~ ., data=train, ntree=1000, mtry=123, importance=TRUE)
mypredict <- predict(rf.fit,test)
table(observed=test$outcome,predicted=mypredict)
confusionMatrix(mypredict , test$outcome) 
acc=length(which(mypredict  == test$outcome))/length(test$outcome)
acc
#=========================================================================================
#gradient boosting
#=========================================================================================
train$outcome = as.numeric(train$outcome)-1
set.seed(17)
gbm.fit <- gbm(formula = outcome ~ .,distribution = "bernoulli",data = train, n.trees = 1000, 
               interaction.depth = 5, shrinkage = 0.001,cv.folds = 10) 
myprediction<-predict(gbm.fit,test,n.trees = 1000,type = "response")
mypred = as.factor(ifelse(myprediction>0.5,1,0))
confusionMatrix(mypred, test$outcome) 
##########################################################################################
#=========================================================================================
#2019
#=========================================================================================
##########################################################################################
#data containing tips before 28/2/2017
data = read.csv("2019/mycurrentdata2019_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
colnames(data)
length(which(data$outcome==1))/nrow(data)
Aux_data = read.csv("Aux_dataNA2019.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
verypast_ind = which(as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) <="2009-1-1")
verypast_name = Aux_data[verypast_ind ,1]
ii = match(verypast_name,data$tip)
length(which(!is.na(ii)))
verypast_ii = ii[which(!is.na(ii))]
length(verypast_ii)
#data$outcome[verypast_ii]
length(which(data$outcome[verypast_ii]==1))
rm_verypast_ii = which(data$outcome[verypast_ii]==0)
verypast_ii = verypast_ii[rm_verypast_ii]
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
alloutliers = c(11048 , 7148 , 7149 ,10551 , 6523 , 6525,  6524 , 6522 , 6521 , 6686, 11169 ,11225 , 5000 , 4725 , 4754 , 4755 , 4807 , 7172 , 7170 , 7173,  7171
                ,6701 , 6707 , 6712,  6714 , 6713 , 6708 , 6709  ,6857 , 7326,  5085 , 6860 , 3234,  3226, 3227 , 3230 , 3231,  3232 , 3229 , 3233,  3228, 15885
                ,15884 , 3742 ,    1 ,10001 ,10002 , 4707, 11951 ,11952)

outliers = alloutliers[1:25]
print(outliers)
#data[outliers,]
data=data[-outliers,]
data_tip = data_tip[-outliers]
length(data_tip)
which(is.na(data),arr.ind = TRUE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load("2019/flutreeHA2019-2.Rdata")  
test_subtree = extract.clade(tree,69893)
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
#ind = which(data$outcome==1)
#data$outcome[ind] = 0
#data$outcome[-ind] =1
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
nrow(test)/(nrow(test)+nrow(train))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#classification with SVM and polynomial kernel
#acc = 0.87 AUC = 0.91
#svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.03125,cost =128,degree=3,coef0=0)
svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.05,cost =64,degree=3,coef0=4)
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
plot(svmmodel.performance,type="l", lwd=1,col="blue", ylab="",xlab="")
par(new=TRUE)
#=========================================================================================
#Random Forest
set.seed(17)
library(randomForest)
# Perform training:
rf.fit = randomForest(outcome ~ ., data=train, ntree=1000, mtry=123, importance=TRUE)
mypredict <- predict(rf.fit,test)
table(observed=test$outcome,predicted=mypredict)
confusionMatrix(mypredict , test$outcome) 
acc=length(which(mypredict  == test$outcome))/length(test$outcome)
acc
#=========================================================================================
#gradient boosting
#=========================================================================================
train$outcome = as.numeric(train$outcome)-1
set.seed(17)
gbm.fit <- gbm(formula = outcome ~ .,distribution = "bernoulli",data = train, n.trees = 1000, 
               interaction.depth = 5, shrinkage = 0.001,cv.folds = 10) 
myprediction<-predict(gbm.fit,test,n.trees = 1000,type = "response")
mypred = as.factor(ifelse(myprediction>0.5,1,0))
confusionMatrix(mypred, test$outcome) 
##########################################################################################
#=========================================================================================
#2018
#=========================================================================================
#source("/Users/maryam/Desktop/Research/NewFlu/code/Tree_Statistics.R")
#source("/Users/maryam/Desktop/Research/NewFlu/code/general_functions.R")
setwd("/Users/maryam/Desktop/Research/NewFlu/GISaid/2018/NA")
#data containing tips before 28/2/2017
data = read.csv("2018/mycurrentdata2018_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
colnames(data)
length(which(data$outcome==1))/nrow(data)
Aux_data = read.csv("Aux_dataNA2018.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
verypast_ind = which(as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) <="2008-1-1")
verypast_name = Aux_data[verypast_ind ,1]
ii = match(verypast_name,data$tip)
length(which(!is.na(ii)))
verypast_ii = ii[which(!is.na(ii))]
length(verypast_ii)
#data$outcome[verypast_ii]
length(which(data$outcome[verypast_ii]==1))
rm_verypast_ii = which(data$outcome[verypast_ii]==0)
verypast_ii = verypast_ii[rm_verypast_ii]
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
#2018
alloutliers = c(2706 , 14436,14438,2994,2092, 2435, 5977,3672, 3664,1836, 817,1726,4653, 9081, 5788 , 5794 , 5789 , 5796, 10268, 13880, 14386,
                13867, 13866 ,13845 , 9828 , 9827 , 4801 , 7195 , 2643 , 1211,  4752,  3665,  5565,  5566 , 3666 , 4604,  8321 , 8184 , 8183,  8167 , 8186,  8185,
                8181 , 8182 , 5064 ,13505 , 1039, 12218 , 3030 , 3211)

outliers = alloutliers[1:25]
print(outliers)
#data[outliers,]
data=data[-outliers,]
data_tip = data_tip[-outliers]
length(data_tip)
which(is.na(data),arr.ind = TRUE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load("2018/flutreeHA2018-2.Rdata")  
test_subtree = extract.clade(tree,60528)
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
#ind = which(data$outcome==1)
#data$outcome[ind] = 0
#data$outcome[-ind] =1
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
nrow(test)/(nrow(test)+nrow(train))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#classification with SVM and polynomial kernel
svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.1,cost =128,degree=3,coef0=1)
print(svm.fit)
mypred=predict(svm.fit,test)
table(mypred, test$outcome) 
sum(mypred==test$outcome)/nrow(test)
confusionMatrix(mypred, test$outcome) 
#=========================================================================================
#Random Forest
set.seed(17)
library(randomForest)
# Perform training:
rf.fit = randomForest(outcome ~ ., data=train, ntree=500, importance=TRUE)
mypredict <- predict(rf.fit,test)
table(observed=test$outcome,predicted=mypredict)
confusionMatrix(mypredict , test$outcome) 
acc=length(which(mypredict  == test$outcome))/length(test$outcome)
acc
#=========================================================================================
#gradient boosting
#=========================================================================================
train$outcome = as.numeric(train$outcome)-1
set.seed(17)
gbm.fit <- gbm(formula = outcome ~ .,distribution = "bernoulli",data = train, n.trees = 500, 
               interaction.depth = 3, shrinkage = 0.1,cv.folds = 5) 
myprediction<-predict(gbm.fit,test,n.trees = 500,type = "response")
mypred = as.factor(ifelse(myprediction>0.5,1,0))
confusionMatrix(mypred, test$outcome) 
##########################################################################################
#=========================================================================================
#2017
#=========================================================================================
#data containing tips before 28/2/2017
data = read.csv("2017/mycurrentdata2017_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
colnames(data)
length(which(data$outcome==1))/nrow(data)
Aux_data = read.csv("Aux_dataNA2017.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
verypast_ind = which(as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) <="2007-1-1")
verypast_name = Aux_data[verypast_ind ,1]
ii = match(verypast_name,data$tip)
length(which(!is.na(ii)))
verypast_ii = ii[which(!is.na(ii))]
length(verypast_ii)
#data$outcome[verypast_ii]
length(which(data$outcome[verypast_ii]==1))
rm_verypast_ii = which(data$outcome[verypast_ii]==0)
verypast_ii = verypast_ii[rm_verypast_ii]
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
alloutliers = c(917,  916,  918,  915, 4531, 2077 , 491,  482 , 483 ,1821 ,5225, 5224, 2766, 2559 ,2439, 3110 ,2611, 2526,
                4665, 4666, 5038 ,3986 ,5332 , 213,  182, 5041, 5124 ,3850, 4713 ,4716 ,4717 ,4714, 4715, 4490 ,5359, 1298,
                1299, 3795,  504 , 500, 4112, 2145 ,4508, 4495 , 190, 2643, 2540 ,2472, 2542, 4476)

outliers = alloutliers[1:25]
print(outliers)
#data[outliers,]
data=data[-outliers,]
data_tip = data_tip[-outliers]
length(data_tip)
which(is.na(data),arr.ind = TRUE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load("2017/flutreeHA2017-2.Rdata")  
test_subtree = extract.clade(tree,47472)
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
#ind = which(data$outcome==1)
#data$outcome[ind] = 0
#data$outcome[-ind] =1
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
nrow(test)/(nrow(test)+nrow(train))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.005,cost =256,degree=4,coef0=0.5)
print(svm.fit)
mypred=predict(svm.fit,test)
table(mypred, test$outcome) 
sum(mypred==test$outcome)/nrow(test)
confusionMatrix(mypred, test$outcome) 
svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
svmmodel.probs<-attr(svmmodel.predict,"decision.values")
svmmodel.class<-predict(svm.fit,test,type="class")
svmmodel.outcome<-test$outcome
#=========================================================================================
#Random Forest
set.seed(17)
library(randomForest)
# Perform training:
rf.fit = randomForest(outcome ~ ., data=train, ntree=1000, mtry=50, importance=TRUE)
mypredict <- predict(rf.fit,test)
table(observed=test$outcome,predicted=mypredict)
confusionMatrix(mypredict , test$outcome) 
acc=length(which(mypredict  == test$outcome))/length(test$outcome)
acc
#=========================================================================================
#gradient boosting
#=========================================================================================
train$outcome = as.numeric(train$outcome)-1
set.seed(17)
gbm.fit <- gbm(formula = outcome ~ .,distribution = "bernoulli",data = train, n.trees = 1000, 
               interaction.depth = 3, shrinkage = 0.1,cv.folds = 10) 
myprediction<-predict(gbm.fit,test,n.trees = 1000,type = "response")
mypred = as.factor(ifelse(myprediction>0.5,1,0))
confusionMatrix(mypred, test$outcome) 
##########################################################################################
#=========================================================================================
#2016  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#data containing tips before 28/2/2017
data = read.csv("2016/mycurrentdata2016_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
colnames(data)
length(which(data$outcome==1))/nrow(data)
Aux_data = read.csv("Aux_dataNA2016.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
verypast_ind = which(as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) <="2006-1-1")
verypast_name = Aux_data[verypast_ind ,1]
ii = match(verypast_name,data$tip)
length(which(!is.na(ii)))
verypast_ii = ii[which(!is.na(ii))]
length(verypast_ii)
#data$outcome[verypast_ii]
length(which(data$outcome[verypast_ii]==1))
rm_verypast_ii = which(data$outcome[verypast_ii]==0)
verypast_ii = verypast_ii[rm_verypast_ii]
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
#2016
alloutliers = c(1233 , 197 , 163,  414 ,1007 ,1014 ,1003, 1075, 1067, 1073, 1058, 2123 ,2136, 2138, 2124, 2126 ,2122, 1416 , 973, 5453, 4973 
                ,2072, 4368 ,4372 ,4373 ,5029, 4972, 4364 ,4362 ,4361, 4889 ,3793, 5700, 5721 ,5886, 2750, 4439 ,3518 ,1936 ,5321, 4383 ,4385, 4384
                ,4319, 4458, 4477 ,4447 ,4450, 6439, 6467)

outliers = alloutliers[1:25]
print(outliers)
#data[outliers,]
data=data[-outliers,]
data_tip = data_tip[-outliers]
length(data_tip)
which(is.na(data),arr.ind = TRUE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load("2016/flutreeHA2016-2.Rdata")  
test_subtree = extract.clade(tree,30199)
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
nrow(test)/(nrow(test)+nrow(train))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#classification
svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.03125,cost =128,degree=3,coef0=0)
print(svm.fit)
mypred=predict(svm.fit,test)
table(mypred, test$outcome) 
sum(mypred==test$outcome)/nrow(test)
confusionMatrix(mypred, test$outcome) 
svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
svmmodel.probs<-attr(svmmodel.predict,"decision.values")
svmmodel.class<-predict(svm.fit,test,type="class")
svmmodel.outcome<-test$outcome
#=========================================================================================
#Random Forest
set.seed(17)
library(randomForest)
# Perform training:
rf.fit = randomForest(outcome ~ ., data=train, ntree=500, mtry=20, importance=TRUE)
mypredict <- predict(rf.fit,test)
table(observed=test$outcome,predicted=mypredict)
confusionMatrix(mypredict , test$outcome) 
acc=length(which(mypredict  == test$outcome))/length(test$outcome)
acc
#=========================================================================================
#gradient boosting
#=========================================================================================
train$outcome = as.numeric(train$outcome)-1
set.seed(17)
gbm.fit <- gbm(formula = outcome ~ .,distribution = "bernoulli",data = train, n.trees = 1000, 
               interaction.depth = 3, shrinkage = 0.1,cv.folds = 10) 
myprediction<-predict(gbm.fit,test,n.trees = 1000,type = "response")
mypred = as.factor(ifelse(myprediction>0.5,1,0))
confusionMatrix(mypred, test$outcome) 
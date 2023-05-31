#predicting the successful strains in 2020
#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")
#BiocManager::install("Biostrings")
#BiocManager::install("GenomicAlignments")

library(e1071)
library(ape)
library(phangorn)
library(apTreeshape)
library(DMwR)
library(outliers)
library(ggplot2)
library("ROCR")
library(caret)
library("gbm")
library("Biostrings")
source("~/code/general_functions.R")
source("~/code/Tree_Statistics.R")
source("~/code/VaccineDist.R")
##IMPORTANT##########################
###Make sure to switch the lables of the model if you swapped the labels during the training!
###########################################################
#train the model
data = read.csv("~/mycurrentdata2020_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
length(which(data$outcome==1))/nrow(data)
Aux_data = read.csv("~/Aux_dataNA2020.csv")
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
load("~/flutree2020-2.Rdata")
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
###########################################################
data_p = read.csv("~/myfutdata_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data_p = data_p[,2:ncol(data_p)]
length(which(data_p$outcome==1))/nrow(data_p)
fut_tip = data_p$tip
data_p = data_p[,-80]
data_p = scaleData(data_p)
data_predict = data_p
#the best model
SVM_pred=predict(svm.fit,data_predict)
successful_tips_ind  = which(SVM_pred == 1)
length(successful_tips_ind )
successful_tips = fut_tip[successful_tips_ind]
#read the epitope distance of the tips
epi_Dis = read.csv("~/epitopeDis.csv",sep= ",",header=T,stringsAsFactors=FALSE)
epi_Dis = epi_Dis[,2:ncol(epi_Dis)]
tip_ind = match(successful_tips,epi_Dis[,3])
length(which(is.na(tip_ind)))
epi_Dis_tips = epi_Dis[tip_ind ,]
plot(density(epi_Dis_tips[,2]))
density(epi_Dis_tips[,2])
median(epi_Dis_tips[,2])
#===================================================================================
#===================================================================================
#we choose our vaccine among the set of sequences after 2019/1
L=sapply(epi_Dis_tips$tip_name , function(x) strsplit(x,"--"))
head(L)
Date = sapply(L, function(x) x[length(x)])
dd = as.Date(Date,"%m/%d/%Y")
ind = which(dd>="2019-1-1")

fut_table = epi_Dis_tips[ind,]
plot(density(fut_table[,2]))
density(fut_table[,2])
median(fut_table[,2])

vaccine_tip = which(fut_table[,2] >= 6.598476 -0.01 & fut_table[,2] <= 6.598476 +0.01)
vaccine_tip_name = fut_table[vaccine_tip,3]
vaccine_tip_name
#===================================================================================
dis_2020=numeric()
for(mytip in vaccine_tip_name){
  tipinf= getVaccineDist_new(mytip, Aux_data,hdata, Pdata,Pdata_2020, year = 2020)
  dis_2020 = c(dis_2020,tipinf$tipinfo)
}
dis_2020
boxplot(dis_2020)
#===================================================================================
#vaccine suggested for 2020-2021: A/Hong Kong/45/2019 (H3N2)-like virus
#compute the distance of our suggested candidates to WHO vaccine
succvaccineDist=getseqDist(vaccine_tip_name, Aux_data,hdata, Pdata, sequence = "A/HongKong/45/2019--EPI_ISL_410589--A_/_H3N2--12/24/2018")
succvaccineDist

#===================================================================================
#random seqs
set.seed(1717)
ind_random = sample(length(successful_tips), length(vaccine_tip_name),replace = FALSE)
random_tips = successful_tips[ind_random]
dis_random_2020=numeric()
for(mytip in random_tips){
  tipinf= getVaccineDist_new(mytip, Aux_data,hdata, Pdata,Pdata_2020, year = 2020)
  dis_random_2020 = c(dis_random_2020,tipinf$tipinfo)
}
dis_random_2020
boxplot(dis_random_2020)

require(e1071)
require(treeWAS)
require(ape)
require(phytools)
require(ROCR)
require(seqinr)
require(ggplot2)

source("~/code/Tree_Statistics.R")
source("~/code/general_functions.R")

# Load data
countrycodes<-read.csv("~/General_data/Countrycode_tips.csv")
data = read.csv("~/2020/mycurrentdata2020_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
Aux_data = read.csv("~2020/Aux_dataNA2020.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]
load("~/2020/flutree2020-2.Rdata")

# Load prediction data
data_p = read.csv("~/2020/myfutdata_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data_p = data_p[,2:ncol(data_p)]
fut_tip = data_p$tip
data_p = data_p[,-80]
data_p = scaleData(data_p)
data_predict = data_p

## Remove strains before 2010
verypast_name = Aux_data[which(as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) <="2010-1-1") ,1]
ii = match(verypast_name,data$tip)
verypast_ii = ii[which(!is.na(ii))]
data = data[-verypast_ii,]

## Remove tip names from data and remove duplicates
data_tip = data$tip
data = data[,-c(80)]
dup_ind = which(duplicated(data))
data = unique(data)
data_tip = data_tip[-dup_ind]

##remove the outliers
set.seed(17)
alloutliers = c(842,1133,880,1072, 306,14 ,115 ,104 ,108,118 , 102 , 110 ,111,107,13015 ,339,341,485,486, 15022, 15033,
                18204, 18584,5133,4128, 21863, 21124,  8611,  6178 ,18498, 18455, 18595, 18577 ,18583 ,18800 ,21125, 21122 ,21121 ,26015, 21507 ,22478, 17370,
                17369 ,26396, 13019 ,13020, 13038, 18684, 18742 ,21758)
outliers = alloutliers[1:25]
data=data[-outliers,]
data_tip = data_tip[-outliers]

### downsample tree by 50% by random sequences by proportion in each region 
success_res<-c(0,0,0)
countrycodes<-countrycodes[countrycodes$TipLabels_Data %in% data_tip,]
samplingprops<-c(0.5,0.2,0.1)

for (n in 1:length(samplingprops)){
  SeqCounts<-data.frame(Region=names(table(countrycodes$region)),
                        No.seqs=round(as.numeric(table(countrycodes$region)*samplingprops[n])))
  
  final_svm<-data.frame(X=0,Y=0,rep=0,auc=0,acc=0)
  times<-0
  repeat {
    
    seqs.in.tree<-character()
    for (k in 1:nrow(SeqCounts)){
      seqs.in.tree<-c(seqs.in.tree,
                      countrycodes$TipLabels_Data[sample(which(countrycodes$region==SeqCounts$Region[k]),
                                                         SeqCounts$No.seqs[k])])
    }
    seqs.in.tree<-c(seqs.in.tree,WHOstrains)
    seqs.in.tree<-unique(seqs.in.tree)
    newtree<-keep.tip(tree,seqs.in.tree)
    newdata_tip<-data_tip[data_tip %in% seqs.in.tree]
    newdata<-data[which(data_tip %in% seqs.in.tree),]
    
    #divide the data to test and train -> test on a test clade and train on other clades
    tipsintree<-testclade_taxa[which(testclade_taxa %in% newtree$tip.label)]
    test_subtree = extract.clade(newtree,findMRCA(newtree,tipsintree))
    test_names = test_subtree$tip.label
    test_ind = match(test_names,newdata_tip)
    testi = test_ind[which(!is.na(test_ind))]
    ind = which(newdata$outcome==1)
    
    newdata$outcome[ind] = 0
    newdata$outcome[-ind] = 1
    newdata = scaleData(newdata)
    test = newdata[testi,]
    test = test[sample(seq_len(nrow(test)), size = nrow(test)),]
    train = newdata[-testi,]
    train = train[sample(seq_len(nrow(train)), size = nrow(train)),]
    
    ## Run support vector machine (SVM) and show prediction
    svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.005,cost =128,degree=5,coef0=0)
    svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
    svmmodel.probs<-attr(svmmodel.predict,"decision.values")
    svmmodel.class<-predict(svm.fit,test,type="class")
    svmmodel.outcome<-test$outcome
    
    #roc analysis for test data
    svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.outcome)
    svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
    if (svmmodel.performance@y.values[[1]][100]>svmmodel.performance@x.values[[1]][100]){
      times<-times+1
      svm_res<-data.frame(X=svmmodel.performance@x.values,Y=svmmodel.performance@y.values,rep=times,auc=0,acc=0)
      svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
      colnames(svm_res)<-c("X","Y","rep","auc","acc")
      mypred=predict(svm.fit,test)
      svm_res$acc<-round(sum(mypred==test$outcome)/nrow(test),2)
      svm_res$auc<-round(svmmodel.auc,2)
      final_svm<-rbind(final_svm,svm_res)
      
      #the best model
      SVM_pred=predict(svm.fit,data_predict)
      successful_tips = fut_tip[which(SVM_pred == 0)]
      success_res<-rbind(success_res,c(length(which(successful_tips %in% successfulStrains))/length(successful_tips),
                                       samplingprops[n],length(which(WHOstrains %in% successful_tips))))
      
    }
    if (times==100){
      break
    }
  }
  final_svm<-final_svm[-1,]
  
  pdf(file = paste0("ROC_plot_downsampling_",samplingprops[n],".pdf"),width = 7,height=6)
  p<-ggplot(final_svm,aes(x=X,y=Y,group=rep)) + geom_line(color="red",alpha=0.4) + theme_bw() + 
    ylab("True Positive Rate") + xlab("False Positive Rate") + geom_abline(color="grey82")
  p<-p + annotate(geom="text", x=0.8, y=0.1, label=paste0("AUC=",median(final_svm$auc)," SD +-",round(sd(final_svm$auc),2),
                                                          "\nAccuracy=",median(final_svm$acc)," SD +-",round(sd(final_svm$acc),2)))
  plot(p)
  dev.off()
}

### downsample tree by 50% by random sequences by proportion in each region 
countrycodes<-countrycodes[countrycodes$TipLabels_Data %in% data_tip,]
SeqCounts<-read.csv("~/General_data/pop_proportions.csv")
SeqCounts$Proportion<-round(SeqCounts$Proportion*6087)


final_svm<-data.frame(X=0,Y=0,rep=0,auc=0,acc=0)
times<-0
repeat {
  
  seqs.in.tree<-character()
  for (k in 1:nrow(SeqCounts)){
    seqs.in.tree<-c(seqs.in.tree,
                    countrycodes$TipLabels_Data[sample(which(countrycodes$region==SeqCounts$Region[k]),
                                                       SeqCounts$Proportion[k])])
  }
  seqs.in.tree<-c(seqs.in.tree,WHOstrains)
  seqs.in.tree<-unique(seqs.in.tree)
  newtree<-keep.tip(tree,seqs.in.tree)
  newdata_tip<-data_tip[data_tip %in% seqs.in.tree]
  newdata<-data[which(data_tip %in% seqs.in.tree),]
  
  #divide the data to test and train -> test on a test clade and train on other clades
  tipsintree<-testclade_taxa[which(testclade_taxa %in% newtree$tip.label)]
  test_subtree = extract.clade(newtree,findMRCA(newtree,tipsintree))
  test_names = test_subtree$tip.label
  test_ind = match(test_names,newdata_tip)
  testi = test_ind[which(!is.na(test_ind))]
  ind = which(newdata$outcome==1)
  
  newdata$outcome[ind] = 0
  newdata$outcome[-ind] = 1
  newdata = scaleData(newdata)
  test = newdata[testi,]
  test = test[sample(seq_len(nrow(test)), size = nrow(test)),]
  train = newdata[-testi,]
  train = train[sample(seq_len(nrow(train)), size = nrow(train)),]
  
  ## Run support vector machine (SVM) and show prediction
  svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.005,cost =128,degree=5,coef0=0)
  svmmodel.predict<-predict(svm.fit, newdata = test,decision.values=TRUE)
  svmmodel.probs<-attr(svmmodel.predict,"decision.values")
  svmmodel.class<-predict(svm.fit,test,type="class")
  svmmodel.outcome<-test$outcome
  
  #roc analysis for test data
  svmmodel.prediction<-prediction(svmmodel.probs,svmmodel.outcome)
  svmmodel.performance<-performance(svmmodel.prediction,"tpr","fpr")
  if (svmmodel.performance@y.values[[1]][100]>svmmodel.performance@x.values[[1]][100]){
    times<-times+1
    svm_res<-data.frame(X=svmmodel.performance@x.values,Y=svmmodel.performance@y.values,rep=times,auc=0,acc=0)
    svmmodel.auc<-performance(svmmodel.prediction,"auc")@y.values[[1]]
    colnames(svm_res)<-c("X","Y","rep","auc","acc")
    mypred=predict(svm.fit,test)
    svm_res$acc<-round(sum(mypred==test$outcome)/nrow(test),2)
    svm_res$auc<-round(svmmodel.auc,2)
    final_svm<-rbind(final_svm,svm_res)
    
    #the best model
    SVM_pred=predict(svm.fit,data_predict)
    successful_tips = fut_tip[which(SVM_pred == 0)]
    success_res<-rbind(success_res,c(length(which(successful_tips %in% successfulStrains))/length(successful_tips),
                                     1,length(which(WHOstrains %in% successful_tips))))
    
  }
  if (times==100){
    break
  }
}
final_svm<-final_svm[-1,]

pdf(file = "ROC_plot_downsampling_propPop.pdf",width = 7,height=6)
p<-ggplot(final_svm,aes(x=X,y=Y,group=rep)) + geom_line(color="red",alpha=0.4) + theme_bw() + 
  ylab("True Positive Rate") + xlab("False Positive Rate") + geom_abline(color="grey82")
p<-p + annotate(geom="text", x=0.8, y=0.1, label=paste0("AUC=",median(final_svm$auc)," SD +-",round(sd(final_svm$auc),2),
                                                        "\nAccuracy=",median(final_svm$acc)," SD +-",round(sd(final_svm$acc),2)))
plot(p)
dev.off()




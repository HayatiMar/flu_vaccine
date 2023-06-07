require(e1071)
require(treeWAS)
require(ape)
require(phytools)
require(ROCR)
require(seqinr)
require(ggplot2)
source("code/Tree_Statistics.R")
source("code/general_functions.R")
source("code/msvmRFE.R")


# Load data
countrycodes<-read.csv("General_data/Countrycode_tips.csv")
SeqCounts<-read.csv("General_data/pop_proportions.csv")
data = read.csv("2020/mycurrentdata2020_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
Aux_data = read.csv("2020/Aux_dataNA2020.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]
load("2020/flutree2020-2.Rdata")

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
countrycodes<-countrycodes[countrycodes$TipLabels_Data %in% data_tip,]
SeqCounts$Proportion<-round(SeqCounts$Proportion*6087)

times<-0
feature_order<-matrix(0,ncol=79,nrow=100)
colnames(feature_order)<-colnames(data)[1:79]
repeat {
  
  times<-times+1
  seqs.in.tree<-character()
  for (k in 1:nrow(SeqCounts)){
    seqs.in.tree<-c(seqs.in.tree,
                    countrycodes$TipLabels_Data[sample(which(countrycodes$region==SeqCounts$Region[k]),
                                                       SeqCounts$Proportion[k])])
  }
  seqs.in.tree<-unique(seqs.in.tree)
  newtree<-keep.tip(tree,seqs.in.tree)
  newdata_tip<-data_tip[data_tip %in% seqs.in.tree]
  newdata<-data[which(data_tip %in% seqs.in.tree),]
  
  #divide the data to test and train -> test on a test clade and train on other clades
  tipsintree<-testclade_taxa[which(testclade_taxa %in% newtree$tip.label)]
  test_subtree = extract.clade(newtree,76876)
  test_names = test_subtree$tip.label
  test_ind = match(test_names,newdata_tip)
  testi = test_ind[which(!is.na(test_ind))]
  ind = which(newdata$outcome==1)
  
  newdata$outcome[ind] = 0
  newdata$outcome[-ind] = 1
  newdata = scaleData(newdata)
  train = newdata[-testi,]
  train = train[sample(seq_len(nrow(train)), size = nrow(train)),]
  train<-cbind(train$outcome,train[,1:79])
  
  feature_selection<-svmRFE(train, k=10)
  feature_order[times,]<-feature_selection
  
  if (times==100){
    break
  }
}

write.csv(feature_order,"feature_order_100repeats_downsample.csv",row.names = F)
feature_order<-feature_order[,order(match(colnames(feature_order),feature_all$feature))]
reshape2::melt(feature_order,'col','row')
melted_feature<-reshape2::melt(as.matrix(feature_order))
melted_feature$Var2<-factor(melted_feature$Var2,levels = unique(melted_feature$Var2))

pdf("FeatureSelection100downsamples.pdf",width = 15, height = 6)
ggplot(melted_feature, aes(x = Var2, y = value)) + geom_boxplot() + ylab("Feature selection (N = 79)") +
  xlab ("Feature") + theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
dev.off()



#library(devtools)
#install_git("http://github.com/bdearlove/treeImbalance.git")
library(treeImbalance)
source("~/Github/Tree_Statistics.R")
library(phyloTop)
#install.packages("remotes")
#remotes::install_github("Leonardini/treeCentrality")
library(treeCentrality)
library(corrplot)

#compute LBI function
#This calculates the local branching index (LBI), defined by Neher, Russell and Shraiman (2014) 
#as the the tree length surrounding a given node or tip exponentially discounted with increasing distance 
#from that node.

Featurs=function(tr){
  features=numeric(40)
  dd = get_pairwise_distances(tr,c(1:length(tr$tip.label)),c(length(tr$tip.label):1))
  Ntip=length(tr$tip.label)
  features[1]=Ntip
  features[2]=sackin(as.treeshape(tr),"pda")
  features[3]=colless(as.treeshape(tr),"pda")
  features[4]=var(node.depth.edgelength(tr)[1:length(tr$tip.label)])
  features[5]=computeI2(tr)
  features[6]=computeB1(tr)
  features[7]=computeB2(tr)
  features[8]=avgLadder(tr, normalise = TRUE)
  features[9]=ILnumber(tr, normalise = TRUE)
  features[10]=pitchforks(tr, normalise = TRUE)
  features[11]=maxHeight(tr, normalise = TRUE)
  features[12]=computeMaxWidth(tr)
  features[13]=computeDelW(tr)
  features[14]=computeStairs1(tr)
  features[15]=computeStairs2(tr)
  features[16]=computeCherries(tr, DOUBLE = FALSE)
  features[17]=computeCherries(tr, DOUBLE = TRUE)
  features[18]=BS(tr)
  features[19]=descinm(tr,m=2)$W
  features[20]=getstattest(tr)$W
  features[21]=skewness(i_bl(tr))
  features[22]=kurtosis(i_bl(tr))
  features[23]=tips_pairwise_distance(tr)
  features[24]=tips_pairwise_distance_max(tr)
  features[25:29]=computeNetworkStats(tr, weight = FALSE, meanpath = FALSE, maxOnly = TRUE)
  features[30:34]=computeNetworkStats(tr, weight = TRUE, meanpath = FALSE, maxOnly = TRUE)
  features[35:38]=computeSpectralStats(tr, weight = c(FALSE, TRUE), adj = c(FALSE,TRUE), norm = FALSE, dist = FALSE, full = FALSE,
                                       maxOnly = TRUE, unitMean = FALSE)
  features[39] = diversificationRate(tr)
  features[40] = mean(lbi(tr,tau = 0.0625*mean(dd)))
  return(features)
}
#random trees
set.seed(1717)
df_features = numeric()
for (i in c(1:200)) {
  x =  sample(c(5:100),1)
  tr = rtree(x)
  df_features = rbind(df_features,Featurs(tr))
}
colnames(df_features) <- c("numberTips","sackin",
                           "colless","Variance","I2","B1","B2","avgLadder","ILnumber","pitchforks",
                           "maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries","DoubleCherries","BS","descinm",
                           "getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist","diameter", "WienerIndex", 
                           "betweenness", "closeness", "eigenvector","diameterW", "WienerIndexW", "betweennessW", "closenessW", 
                           "eigenvectorW","minAdj","maxAdj","minLap","maxLap","DivR","LBI")

dd = cor(df_features)
corrplot(dd, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,tl.cex = 0.6)
#=======================================================================================
#HA tree
load("~/2020/flutree2020-2.Rdata")
set.seed(1234)
df_features = numeric()
for (i in c(1:200)) {
  x =  sample((length(tree$tip.label)+1):(2*length(tree$tip.label)-1),1)
  print(x)
  tr = extract.clade(tree,x)
  if(length(tr$tip.label) >=5 & length(tr$tip.label) <= 200){
  df_features = rbind(df_features,Featurs(tr))
  }
}
colnames(df_features) <- c("numberTips","sackin",
                           "colless","Variance","I2","B1","B2","avgLadder","ILnumber","pitchforks",
                           "maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries","DoubleCherries","BS","descinm",
                           "getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist","diameter", "WienerIndex", 
                           "betweenness", "closeness", "eigenvector","diameterW", "WienerIndexW", "betweennessW", "closenessW", 
                           "eigenvectorW","minAdj","maxAdj","minLap","maxLap","DivR","LBI")

dd = cor(df_features)
corrplot(dd, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,tl.cex = 0.6)
#=======================================================================================
#NA tree
load("~/2020flutree_NA_2020-2.Rdata")
set.seed(1234)
df_features = numeric()
for (i in c(1:200)) {
  x =  sample((length(tree$tip.label)+1):(2*length(tree$tip.label)-1),1)
  print(x)
  tr = extract.clade(tree,x)
  if(length(tr$tip.label) >=5 & length(tr$tip.label) <= 200){
    df_features = rbind(df_features,Featurs(tr))
  }
}
colnames(df_features) <- c("numberTips_NA","sackin_NA",
                           "colless_NA","Variance_NA","I2_NA","B1_NA","B2_NA","avgLadder_NA","ILnumber_NA","pitchforks_NA",
                           "maxHeight_NA","MaxWidth_NA","DelW_NA","Stairs1_NA","Stairs2_NA","Cherries_NA","DoubleCherrie_NAs","BS_NA","descinm_NA",
                           "getstattest_NA","skewness_NA","kurtosis_NA","MeanPairwiseDist_NA","MaxPairwiseDist_NA","diameter_NA", "WienerIndex_NA", 
                           "betweenness_NA", "closeness_NA", "eigenvector_NA","diameterW_NA", "WienerIndexW_NA", "betweennessW_NA", "closenessW_NA", 
                           "eigenvectorW_NA","minAdj_NA","maxAdj_NA","minLap_NA","maxLap_NA","DivR_NA","LBI_NA")

dd = cor(df_features)
corrplot(dd, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,tl.cex = 0.6)

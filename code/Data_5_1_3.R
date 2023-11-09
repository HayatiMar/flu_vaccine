library(ape)
library(phangorn)
library(apTreeshape)
library(igraph)
library(phyloTop)
library(moments)
library(devtools)
devtools::install_github('Leonardini/treeCentrality')
library(treeCentrality)
library(seqinr)
source("code/Tree_Statistics.R")
source("code/general_functions.R")

load("2020/flutreeHA2020-2.Rdata")
Aux_data = read.csv("2020/Aux_dataHA2020.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]
treeinds=match(tree$tip.label, Aux_data[,1])
Aux_data = Aux_data[treeinds,]
all(tree$tip.label == Aux_data[treeinds,1])


myclades=getClades2(tree,MinTotalSize=5, MinTrimSize=5, OverlapCutoff=1, TimeFrame=3)
Final_trimmedClades=myclades$trimclades[myclades$rejected==0]
Final_trimmedClades_root=as.numeric(names(which(myclades$rejected==0)))
length(Final_trimmedClades_root)
allHeights=node.depth.edgelength(tree); max(allHeights)  #
hdata=data.frame(tiplab=tree$tip.label, height=allHeights[1:length(tree$tip.label)])
Labels_tab = GetLabels(tree)
length(which(Labels_tab[,"allTrimmedClades_3_size"]>=5))
#remove the very recent clade (we use them in prediction task)
res=numeric()
for(i in 1:length(Final_trimmedClades_root)){
  #print(i)
  if( max(allHeights)-allHeights[Final_trimmedClades_root[i]] <=4 )
    {res=c(res,i)}
}
length(res)
Final_root2020=Final_trimmedClades_root[-res] 
length(Final_root2020)
#remove the very big clades
ind = match(Final_root2020, Labels_tab[,1])
ind_1 = which(Labels_tab[ind,"allTrimmedClades_4_size"]<=875)  
length(ind_1)
ind_2 = match(Labels_tab[ind,][ind_1,1],Final_root2020)
Final_root2020 = Final_root2020[ind_2]
length(Final_root2020)

#=================extract the features=================================================
Data_Clade=matrix(0,length(Final_trimmedClades),41)

for(i in match(Final_root2020,Final_trimmedClades_root)){
  print(Final_trimmedClades_root[i])
  tr1=extract.clade(tree,Final_trimmedClades_root[i])
  tr=drop.tip(tr1,setdiff(tr1$tip.label,hdata$tiplab[Final_trimmedClades[[i]]]), trim.internal = TRUE)
  root=Final_trimmedClades_root[i]
  Data_Clade[i,]=Clade_featurs(tree,tr1,tr,root)
}
Data_Clade = as.matrix(Data_Clade)
Data_Clade = Data_Clade[-which(Data_Clade[,1]==0),]
#====================================================================================

colnames(Data_Clade) = c("Clade","numberTipsClade","numberTipsTrimmed","sackin",
                 "colless","Variance","I2","B1","B2","avgLadder","ILnumber","pitchforks",
                 "maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries","DoubleCherries","BS","descinm",
                 "getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist","diameter", "WienerIndex", 
                 "betweenness", "closeness", "eigenvector","diameterW", "WienerIndexW", "betweennessW", "closenessW", 
                 "eigenvectorW","minAdj","maxAdj","minLap","maxLap","DivR")

#====================================================================================
#find the labels
allHeights=node.depth.edgelength(tree)
allD=allDescendants(tree) 
TimeFrame=4
nodeids=Data_Clade[,1]
nTips=length(tree$tip.label)
# # tips within the TimeFrame for each node
allTrimmedClades = sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
myTipTimes=allHeights[myTipDes] # here, would need something like allDates[myTipDesc]
return(myTipDes[myTipTimes <= allHeights[x]+TimeFrame]) })
# # sizes of trimmed clades 
allTrimmedSizes_4 = sapply(allTrimmedClades, function(x) length(allTrimmedClades[x]))
Labels=allTrimmedSizes_4/Data_Clade[,3]
Data_Clade = cbind(Data_Clade,cbind(allTrimmedSizes_4,Labels))
df_names = c("Clade","numberTipsClade","numberTipsTrimmed","sackin",
             "colless","Variance","I2","B1","B2","avgLadder","ILnumber","pitchforks",
             "maxHeight","MaxWidth","DelW","Stairs1","Stairs2","Cherries","DoubleCherries","BS","descinm",
             "getstattest","skewness","kurtosis","MeanPairwiseDist","MaxPairwiseDist","diameter", "WienerIndex", 
             "betweenness", "closeness", "eigenvector","diameterW", "WienerIndexW", "betweennessW", "closenessW", 
             "eigenvectorW","minAdj","maxAdj","minLap","maxLap","DivR","numberTipsTrimmed_4","Labels")
names(Data_Clade) = df_names

write.csv(Data_Clade,"Data5_1_3.csv")


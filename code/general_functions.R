
#' @param tree phylogenetic tree in phylo format
#' @param weighted true/false - if true, use the tree's branch lengths. if false, use the number of edges 
#' @return a list, each entry corresponding to a tip. each element of the list is a vector of distances from the tip 
#' each node on the path between the tip and the root, measured in units of branch length by default 
#' @examples 
#' getdtk(rtree(20))
getdtk = function(tree,weighted=TRUE) {
  num.tips=length(tree$tip.label)
  if (weighted==FALSE) tree$edge.length=1+0*tree$edge.length 
  np=nodepath(tree) # list of vectors listing paths from root to tip 1, root to tip 2 etc   
  allHeights=node.depth.edgelength(tree);
  dtk=np
  for (tip in 1:num.tips) {
    dtk[[tip]] = allHeights[tip] - allHeights[ np[[tip]] ]
    names(dtk[[tip]]) = np[[tip]] # names will now be the node ids on the paths from tips to root
  }
  return(dtk)
}


getTipFeatures = function(tree, dtkweighted=TRUE, rescale=TRUE,sig=5*median(tree$edge.length),sizescale=FALSE,df,TimeFrame =3) {
  num.tips=length(tree$tip.label)
  clade_names = df$Clade 
  df$Clade <- NULL# not needed here - here we only want the features, not the node ids
  df$Labels <-NULL
  if (rescale==TRUE) df=scale(df) 
  if (rescale==FALSE) df = as.matrix(df)
  dtk=getdtk(tree,weighted=dtkweighted)
  Heights=node.depth.edgelength(tree)
  # preallocate
  tipdata = matrix(0, nrow=num.tips, ncol=ncol(df))
  for (tip in 1:num.tips) {
    print(tip)
    ind = which(!is.na(match(as.numeric(names(dtk[[tip]])),clade_names)))
    clade_names_ind = match(as.numeric(names(dtk[[tip]])),clade_names)[rev(ind)]
    tip_h = Heights[tip]
    good_clades = which(tip_h - Heights[clade_names[clade_names_ind]] >= TimeFrame)
    if (length( good_clades) >= 1 ){
      ii = clade_names_ind[good_clades]
      MyMatrix=df[ii,]  
      if(length( good_clades) == 1){tipdata[tip,]=MyMatrix}
      if(length( good_clades) > 1){
      arr_i = match(clade_names[ii],names(dtk[[tip]]))
      arr = dtk[[tip]][arr_i]
      arr_t = arr/sig
      coef = trans_interval(arr_t)
      myVec=exp(-(arr_t*coef[1]+coef[2]))
      tipdata[tip,]=myVec %*% MyMatrix
      }
    }
  }
  tipdata=as.data.frame(tipdata)
  colnames(tipdata)=colnames(df)
  tipdata$tip=tree$tip.label
  return(tipdata)
}


#' @param tree phylogenetic tree in phylo format, with data from both past and 'present'
#' @param dtkweighted  Logical: use branch lengths in computing the distances from tip to nodes on path to root?
#' @param sig scale in exponential for distances along path to root
#' @param df dataset
#' @param TimeFrame the time frame that the subtrees are trimmed 
#' @return vector of response variables in the given mode, corresponding to tips from the past only. present tips do not get a response
getResponse = function(tree,dtkweighted=TRUE,sig=5*mean(tree$edge.length), df,TimeFrame =3) { 
  num.tips=length(tree$tip.label)
  clade_names =df$Clade 
  Heights=node.depth.edgelength(tree)
  Heights_tips=node.depth.edgelength(tree)[1: num.tips]
  dtk=getdtk(tree,weighted=dtkweighted) # we don't need all the tips 
  tipresponse=0*(1:num.tips)
  names(tipresponse)=tree$tip.label
  fraction=df$Labels 
  tipresponses=data.frame(tip=tree$tip.label,origheight=Heights_tips, frac=0*(1:length(num.tips)))
  # for each tip we want a weighted sum of the above values on paths from the tip to the root
  for (tip in 1:num.tips) {
    #print(tip)
    ind = which(!is.na(match(as.numeric(names(dtk[[tip]])),clade_names)))
    clade_names_ind = match(as.numeric(names(dtk[[tip]])),clade_names)[rev(ind)]
    tip_h = Heights[tip]
    good_clades = which(tip_h - Heights[clade_names[clade_names_ind ]] >= TimeFrame)
    if (length( good_clades) >= 1 ){
      noderesp.frac = fraction[match(clade_names[clade_names_ind][good_clades],df[,1])]
      ii = clade_names_ind[good_clades]
      arr_i = match(clade_names[ii],names(dtk[[tip]]))
      if(length( good_clades) == 1){tipresponses$frac[tip]=noderesp.frac}
      if(length( good_clades) > 1){
      arr = dtk[[tip]][arr_i]
      arr_t = arr/sig
      coef = trans_interval(arr_t)
      myVec=exp(-(arr_t*coef[1]+coef[2]))
      tipresponses$frac[tip]=myVec %*% noderesp.frac 
      }
    }
  }
  return(tipresponses)
}

#' @param tree is a tree in phylo format
#' @param tr1 is the subtree rooted at node root
#' @param tr is the trimmed subtree rooted at node root
#' @param hdata is a matrix with two columns: the first stores the tip labels and the second column stores the height of each tip
#' @return the features of a given subtree 
#' @examples
Clade_featurs=function(tree,tr1,tr,root,hdata){
  Ntip1=length(tr1$tip.label)
  Ntip2=length(tr$tip.label)
  features=numeric(41)
  features[1]=root
  features[2]=Ntip1
  features[3]=Ntip2
  features[4]=sackin(as.treeshape(tr),"pda")
  features[5]=colless(as.treeshape(tr),"pda")
  features[6]=var(node.depth.edgelength(tr)[1:length(tr$tip.label)])
  features[7]=computeI2(tr)
  features[8]=computeB1(tr)
  features[9]=computeB2(tr)
  features[10]=avgLadder(tr, normalise = TRUE)
  features[11]=ILnumber(tr, normalise = TRUE)
  features[12]=pitchforks(tr, normalise = TRUE)
  features[13]=maxHeight(tr, normalise = TRUE)
  features[14]=computeMaxWidth(tr)
  features[15]=computeDelW(tr)
  features[16]=computeStairs1(tr)
  features[17]=computeStairs2(tr)
  features[18]=computeCherries(tr, DOUBLE = FALSE)
  features[19]=computeCherries(tr, DOUBLE = TRUE)
  features[20]=BS(tr)
  features[21]=descinm(tr,m=2)$W
  features[22]=getstattest(tr)$W
  features[23]=skewness(i_bl(tr))
  features[24]=kurtosis(i_bl(tr))
  features[25]=tips_pairwise_distance(tr)
  features[26]=tips_pairwise_distance_max(tr)
  features[27:31]=computeNetworkStats (tr, weight = FALSE, meanpath = FALSE, maxOnly = TRUE)
  features[32:36]=computeNetworkStats (tr, weight = TRUE, meanpath = FALSE, maxOnly = TRUE)
  features[37:40]=computeSpectralStats(tr, weight = c(FALSE, TRUE), adj = c(FALSE,TRUE), norm = FALSE, dist = FALSE, full = FALSE,
                                       maxOnly = TRUE, unitMean = FALSE)
  features[41] = diversificationRate(tr)
  return(features)
}


#' @param tree is a tree in phylo format
#' @param Final_trimmedClades_root the root of the accepted trimmed clades using getclade function
#' @param Final_trimmedClades the tips include in the accepted trimmed clades using getclade function
#' @param numFeatures the number of features that will be computed 
#' @return a dataframe which each row is the features of each accepted subtrees of a given tree
#' @examples
computeFeatures = function(tree,Final_trimmedClades_root,Final_trimmedClades,numFeatures){
  allHeights=node.depth.edgelength(tree)
  hdata=data.frame(tiplab=tree$tip.label, height=allHeights[1:length(tree$tip.label)])
  Data_Clade=matrix(0,length(Final_trimmedClades_root),numFeatures)
  for(i in 1:length(Final_trimmedClades_root)){
    if(i%%2 == 0 ){print(i)}
    tr1=extract.clade(tree,Final_trimmedClades_root[i])
    tr=drop.tip(tr1,setdiff(tr1$tip.label,hdata$tiplab[Final_trimmedClades[[i]]]), trim.internal = TRUE)
    root=Final_trimmedClades_root[i]
    Data_Clade[i,]=Clade_featurs(tree,tr1,tr,root,hdata)
  }
  Data_Clade = as.matrix(Data_Clade)
  return(Data_Clade)
}


#' @param tree is a tree in phylo format
#' @return a vector giving the size of the subtrees trimmed after different time-frame
#' @examples
GetLabels = function(tree){
  allHeights=node.depth.edgelength(tree)
  allD=allDescendants(tree) 
  nTips=length(tree$tip.label)
  nodeids = (nTips+1):(2*nTips-1)
  # number of tips in the timeframe from the root of the subtree to the (date of the last tip in the tree - 1)
  allTrimmedClades = sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes] 
  return(myTipDes[myTipTimes <= (max(allHeights[myTipDes])-1)])})
  # sizes of trimmed clades 
  allTrimmedSizes = sapply(allTrimmedClades, function(x) length(allTrimmedClades[x]))
  #the size of all the subtrees
  allCladeSizes = sapply(nodeids, function(x) {myTipDesLength=length(allD[[x]][allD[[x]]<=nTips])
  return(myTipDesLength)})
  #the height of the root of the subtrees
  allheightNodes = allHeights[nodeids]
  # the height of last tip of each subtree
  LastTips = sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes] 
  return(max(allHeights[myTipDes]))})
  # number of tips within the TimeFrame for each node
  allTrimmedClades_3= sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes] 
  return(myTipDes[myTipTimes <= allHeights[x]+3]) })
  allTrimmedClades_3_size = sapply(allTrimmedClades_3, function(x) length(allTrimmedClades_3[x])) 
  allTrimmedClades_2= sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes] 
  return(myTipDes[myTipTimes <= allHeights[x]+2]) })
  allTrimmedClades_2_size = sapply(allTrimmedClades_2, function(x) length(allTrimmedClades_2[x]))
  allTrimmedClades_4= sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes]
  return(myTipDes[myTipTimes <= allHeights[x]+4]) })
  allTrimmedClades_4_size = sapply(allTrimmedClades_4, function(x) length(allTrimmedClades_4[x]))
  allTrimmedClades_5= sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes] 
  return(myTipDes[myTipTimes <= allHeights[x]+5]) })
  allTrimmedClades_5_size = sapply(allTrimmedClades_5, function(x) length(allTrimmedClades_5[x]))
  allTrimmedClades_6= sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes]
  return(myTipDes[myTipTimes <= allHeights[x]+6]) })
  allTrimmedClades_6_size = sapply(allTrimmedClades_6, function(x) length(allTrimmedClades_6[x]))
  allTrimmedClades_3.4= sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes] 
  return(myTipDes[myTipTimes <= allHeights[x]+3.4]) })
  allTrimmedClades_3.4_size = sapply(allTrimmedClades_3.4, function(x) length(allTrimmedClades_3.4[x]))
  allTrimmedClades_1.4= sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes]
  return(myTipDes[myTipTimes <= allHeights[x]+1.4]) })
  allTrimmedClades_1.4_size = sapply(allTrimmedClades_1.4, function(x) length(allTrimmedClades_1.4[x]))
  allTrimmedClades_1= sapply(nodeids, function(x) {myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes] 
  return(myTipDes[myTipTimes <= allHeights[x]+1]) })
  allTrimmedClades_1_size = sapply(allTrimmedClades_1, function(x) length(allTrimmedClades_1[x]))
  
  Labels_tab = cbind(nodeids,allheightNodes,LastTips,allCladeSizes, allTrimmedSizes,allTrimmedClades_3_size, allTrimmedClades_2_size,allTrimmedClades_4_size, allTrimmedClades_5_size, allTrimmedClades_6_size, 
                     allTrimmedClades_1_size,allTrimmedClades_1.4_size, allTrimmedClades_3.4_size)
  colnames(Labels_tab) = c("nodeids","allheightNodes","LastTips","allCladeSizes",            
                        "allTrimmedSizes","allTrimmedClades_3_size","allTrimmedClades_2_size","allTrimmedClades_4_size",  "allTrimmedClades_5_size","allTrimmedClades_6_size",
                        "allTrimmedClades_1_size","allTrimmedClades_1.4_size","allTrimmedClades_3.4_size")
  return(Labels_tab)
}

trans_interval = function(arr){
  A = matrix(c(min(arr),max(arr),1,1),nrow=2,ncol=2)
  b = matrix(c(0, 1.5*max(arr)),nrow=2,ncol=1)
  return(solve(A, b))
}


diversificationRate = function(tr){
  T_sum = 0
  for (i in nodepath(tr)) {
    Nedges = length(i)-1
    sum = 0
    for (j in 1:(length(i)-1)) {
      ind = which(apply(tr$edge, 1, function(x) all.equal(x, c(i[j],i[j+1]))) == "TRUE")  
      sum = sum + tr$edge.length[ind]/(2^(j-1))
    }
    sum = sum^-1
    T_sum = T_sum+sum
  }
  T_sum = T_sum/length(nodepath(tr))
  return(T_sum)  
}


scaleData=function(data){
  Ndata=colnames(data)
  scaled.data <- scale(data[,1:(ncol(data)-1)])
  scaled.data=as.data.frame(scaled.data)
  data$outcome=as.factor(data$outcome)
  data=cbind(scaled.data ,data$outcome)
  colnames(data)=Ndata
  return(data)
}



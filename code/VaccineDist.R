library(ape)
library(phangorn)
library(seqinr)
#install.packages("alakazam")
library(alakazam)

load("2020/flutree2020-2.Rdata")  
tree_lab=tree$tip.label
LR=sapply(tree_lab, function(x) strsplit(x,"--"))
R_Lab = sapply(LR, function(x) x[2])
head(R_Lab)

allHeights=node.depth.edgelength(tree)
hdata=data.frame(tiplab=tree$tip.label, height=allHeights[1:length(tree$tip.label)])

Pdata=read.fasta("General_data/Pdata.fasta",seqtype="AA")
PTipLabels_Data = names(Pdata)
LP=sapply(PTipLabels_Data, function(x) strsplit(x,"--"))
P_Lab = sapply(LP, function(x) x[2])
head(P_Lab)

ind=match(R_Lab,P_Lab)
NA_ind=which(is.na(ind))
length(NA_ind)

#tree_lab[NA_ind] 
tree_lab = tree_lab[-NA_ind]

PTipLabels_Data = PTipLabels_Data[ind]
PTipLabels_Data = PTipLabels_Data[-NA_ind]
Pdata = Pdata[ind]
Pdata = Pdata[-NA_ind]
#check
PTipLabels_Data[length(PTipLabels_Data)]
tree_lab[[length(PTipLabels_Data)]]
Pdata[length(PTipLabels_Data)]

Aux_data = read.csv("2020/Aux_dataHA2020.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]

ind = match(tree_lab,Aux_data[,1])
NA_ind=which(is.na(ind))
Aux_data = Aux_data[ind,]

ind = match(tree_lab,hdata[,1])
NA_ind=which(is.na(ind))
hdata = hdata[ind,]


LR=sapply(tree_lab, function(x) strsplit(x,"--"))
R_Lab = sapply(LR, function(x) x[2])
head(R_Lab)

Date = sapply(LR, function(x) x[length(x)])
head(Date)
dd=as.Date(Date,"%m/%d/%Y")
head(dd)
#==============================================================================================
# the HA1 subunit starts at aa 17; this is stated in the NCBI annotation
epitopesites=16+c(2,3,5,25,33,50,53,54,57,62,63,67,75,78,81,82,83,92,94,106,121,122,124,126,131,133,135,
                  137,142,143,144,145,146,155,156,157,158,159,160,163,164,172,173,174,186,188,189,190,192,193,196,197,
                  201,202,207,213,217,222,225,226,227,242,244,248,260,262,271,275,276,278,299,307)



getEpitopeSites <- function(sdata,epitopesites){
  newdata=t(sapply(sdata, function(x) x[epitopesites]))
  newdata=as.AAbin(newdata)
  return(newdata)
}

getVaccineDist <- function(mytip, Aux_data,hdata, Pdata, year = 2019) {
  if (is.character(mytip)) { tipind=match(mytip, hdata$tiplab)} else {tipind=mytip}
  tipTime=hdata$height[tipind]
  tiplabel = hdata$tiplab[tipind]
  tiplab_AA=names(Pdata)[tipind]
  if (year == 2020){low ="2020-2-28" 
  high = "2021-2-28" }
  if (year == 2019){low ="2019-2-28" 
  high = "2020-2-28" }
  if (year == 2018){low ="2018-2-28" 
  high = "2019-2-28" }
  if (year == 2017){low ="2017-2-28" 
  high = "2018-2-28" }
  if (year == 2016){low ="2016-2-28" 
  high = "2017-2-28" }
  reltipindex= which(as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) >= low & as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) <=high)
  if (length(reltipindex)==0) reltipindex=tipind
  # find relevant tips
  #relTips=hdata$tiplab[reltipindex]
  relTips_AA=names(Pdata)[reltipindex]
  #find the AA labels and index of relevant tips
  lmt=length(mytip); lrt=length(relTips_AA)
  # if present, remove everything but the epitope sites
  if (length(Pdata[[1]]) > 72) {Cur_Pdata = getEpitopeSites(Pdata[c(tiplab_AA, relTips_AA)],epitopesites)}
  #find an appropariate format for 
  Cur_Pdata_list = as.list(Cur_Pdata)
  Cur_Pdata_list_new = list()
  for (i in Cur_Pdata_list) {
    nn=character()
    for (j in as.character(i)) {
      nn = paste0(nn,j)
    }
    Cur_Pdata_list_new = c(Cur_Pdata_list_new,nn)
  }                             
  
  mydd = sapply(Cur_Pdata_list_new[2:length(Cur_Pdata_list_new)], function(x) 
    seqDist(Cur_Pdata_list_new[[1]], x, dist_mat = getAAMatrix()))
  #mydd = exp(-mydd/D0)
  # compute relevant distance: for each tip in my clade I want the sum of all its dists to rel tips
  # which is the sum over the row of exp(-Dij/Do)
  tipinfo=sum(mydd)/lrt # should be as many of these as there are tips in my clade; div by lrt to use mean
  # so that I can compare clades with different sizes of past tips
  #names(tipinfo)=mytip
  res = data.frame(mytip, tipinfo)
  rownames(res)=tiplabel
  return(res)
}

Pdata_2020=read.fasta("General_data/Pdata_2020-21.fasta",seqtype="AA")
getVaccineDist_new <- function(mytip, Aux_data,hdata, Pdata,Pdata_2020, year = 2020) {
  if (is.character(mytip)) { tipind=match(mytip, hdata$tiplab)} else {tipind=mytip}
  tipTime=hdata$height[tipind]
  tiplabel = hdata$tiplab[tipind]
  tiplab_AA=names(Pdata)[tipind]
  if (year == 2020){
    #low ="2020-2-28" 
    #high = "2021-2-28" 
    Pdata_new=Pdata_2020
    reltipindex = c(1:length(Pdata_new))
    Pdata = c(Pdata_new,Pdata)
  }
  if (year == 2019){low ="2019-2-28" 
  high = "2020-2-28" }
  if (year == 2018){low ="2018-2-28" 
  high = "2019-2-28" }
  if (year == 2017){low ="2017-2-28" 
  high = "2018-2-28" }
  if (year == 2016){low ="2016-2-28" 
  high = "2017-2-28" }
  #reltipindex= which(as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) >= low & as.Date(as.Date(Aux_data$Date,"%m/%d/%Y")) <=high)
  #if (length(reltipindex)==0) reltipindex=tipind
  # find relevant tips
  #relTips=hdata$tiplab[reltipindex]
  relTips_AA=names(Pdata)[reltipindex]
  #find the AA labels and index of relevant tips
  lmt=length(mytip); lrt=length(relTips_AA)
  # if present, remove everything but the epitope sites
  if (length(Pdata[[1]]) > 72) {Cur_Pdata = getEpitopeSites(Pdata[c(tiplab_AA, relTips_AA)],epitopesites)}
  #find an appropariate format for 
  Cur_Pdata_list = as.list(Cur_Pdata)
  Cur_Pdata_list_new = list()
  for (i in Cur_Pdata_list) {
    nn=character()
    for (j in as.character(i)) {
      nn = paste0(nn,j)
    }
    Cur_Pdata_list_new = c(Cur_Pdata_list_new,nn)
  }                             
  
  mydd = sapply(Cur_Pdata_list_new[2:length(Cur_Pdata_list_new)], function(x) 
    seqDist(Cur_Pdata_list_new[[1]], x, dist_mat = getAAMatrix()))
  #mydd = exp(-mydd/D0)
  # compute relevant distance: for each tip in my clade I want the sum of all its dists to rel tips
  # which is the sum over the row of exp(-Dij/Do)
  tipinfo=sum(mydd)/lrt # should be as many of these as there are tips in my clade; div by lrt to use mean
  # so that I can compare clades with different sizes of past tips
  #names(tipinfo)=mytip
  res = data.frame(mytip, tipinfo)
  rownames(res)=tiplabel
  return(res)
}

getseqDist <- function(mytips, Aux_data,hdata, Pdata, sequence = "A/Hong_Kong/45/2019--EPI_ISL_390308--A_/_H3N2--12/24/2018") {
  if (is.character(mytips)) { tipind=match(mytips, hdata$tiplab)} else {tipind=mytips}
  tipTime=hdata$height[tipind]
  tiplabel = hdata$tiplab[tipind]
  tiplab_AA=names(Pdata)[tipind]

  reltipindex= which(names(Pdata)==sequence)
  #if (length(reltipindex)==0) reltipindex=tipind
  # find relevant tips
  #relTips=hdata$tiplab[reltipindex]
  relTips_AA=names(Pdata)[reltipindex]
  #find the AA labels and index of relevant tips
  lmt=length(mytips); lrt=length(relTips_AA)
  # if present, remove everything but the epitope sites
  if (length(Pdata[[1]]) > 72) {Cur_Pdata = getEpitopeSites(Pdata[c(relTips_AA,tiplab_AA)],epitopesites)}
  #find an appropariate format for 
  Cur_Pdata_list = as.list(Cur_Pdata)
  Cur_Pdata_list_new = list()
  for (i in Cur_Pdata_list) {
    nn=character()
    for (j in as.character(i)) {
      nn = paste0(nn,j)
    }
    Cur_Pdata_list_new = c(Cur_Pdata_list_new,nn)
  }                             
  
  mydd = sapply(Cur_Pdata_list_new[2:length(Cur_Pdata_list_new)], function(x) 
    seqDist(Cur_Pdata_list_new[[1]], x, dist_mat = getAAMatrix()))
  names(mydd) = tiplab_AA
  return(mydd)
}



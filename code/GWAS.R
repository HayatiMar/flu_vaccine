require(e1071)
require(treeWAS)
require(ape)
require(phytools)
require(ROCR)
require(seqinr)
source("~/code/Tree_Statistics.R")
source("~/code/general_functions.R")

# Load data
data = read.csv("~/2020/mycurrentdata2020_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data = data[,2:ncol(data)]
Aux_data = read.csv("~2020/Aux_dataNA2020.csv")
Aux_data = Aux_data[,2:ncol(Aux_data)]
load("~/2020/flutree2020-2.Rdata")

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

#divide the data to test and train -> test on a test clade and train on other clades
test_subtree = extract.clade(tree,76876)
test_names = test_subtree$tip.label
test_ind = match(test_names,data_tip)
testi = test_ind[which(!is.na(test_ind))]
ind = which(data$outcome==1)
data$outcome[ind] = 0
data$outcome[-ind] = 1
data = scaleData(data)
test = data[testi,]
test = test[sample(seq_len(nrow(test)), size = nrow(test)),]
train = data[-testi,]
train = train[sample(seq_len(nrow(train)), size = nrow(train)),]

## Run support vector machine (SVM) and show prediction
svm.fit = svm(data = train, outcome ~ .,kernel ="polynomial",  scale =FALSE,gamma = 0.005,cost =128,degree=5,coef0=0)

# Load future data
data_p = read.csv("~/2020/myfutdata_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data_p = data_p[,2:ncol(data_p)]
fut_tip = data_p$tip
data_p = data_p[,-80]
data_p = scaleData(data_p)
data_predict = data_p

#the best model
SVM_pred=predict(svm.fit,data_predict)
successful_tips = fut_tip[which(SVM_pred == 0)]
nonsuccessful_tips = fut_tip[which(SVM_pred == 1)]
prediction_success<-rbind(cbind(successful_tips,1),cbind(nonsuccessful_tips,0))
colnames(prediction_success)<-c("Name","Success")
write.table(prediction_success,"SuccessfulStrains.txt",row.names = F,quote = F,sep = "\t")

# Run GWAS HA
fastafile<-"Sequences/HA_aligned_clean.fasta"
load("../trees/flutree2020-2.Rdata")
phenofile<-read.table("SuccessfulStrainsNew.txt",header = T)

dna <- read.FASTA(file = fastafile,type = "DNA")
dna<-dna[names(dna) %in% phenofile[,1]]
tree<-keep.tip(tree,names(dna))
mat <- DNAbin2genind(dna)@tab
phenofile <- phenofile[which(phenofile[,1] %in% tree$tip.label),]
phen <- as.vector(as.numeric(unlist(phenofile[,2])))
names(phen) <- phenofile[,1]

out_HA <- treeWAS(snps = mat,
               phen = phen,
               tree = tree,
               seed = 1,test = "terminal",plot.tree = F,plot.null.dist = F)


# Run GWAS NA
fastafile<-"NA_aligned_clean.fasta" ## Download all NA sequences from GISAID
load("~/2020/flutree2020-2.Rdata")
phenofile<-read.table("SuccessfulStrains.txt",header = T)

dna <- read.FASTA(file = fastafile,type = "DNA")
dna<-dna[names(dna) %in% phenofile[,1]]
tree<-keep.tip(tree,names(dna))
mat <- DNAbin2genind(dna)@tab
phenofile <- phenofile[which(phenofile[,1] %in% tree$tip.label),]
phen <- as.vector(as.numeric(unlist(phenofile[,2])))
names(phen) <- phenofile[,1]

out_NA <- treeWAS(snps = mat,
                  phen = phen,
                  tree = tree,
                  seed = 1,test = "terminal",plot.tree = F,plot.null.dist = F)

##significant SNPs

sig.SNPs.NA<-out_NA$treeWAS.combined$treeWAS.combined
sig.SNPs.HA<-out_HA$treeWAS.combined$treeWAS.combined

prop.successful<-length(which(names(which(mat[,which(colnames(mat)=="873.c")]==1)) %in% successful_tips))/length(successful_tips)
prop.nonsuccessful<-length(which(names(which(mat[,which(colnames(mat)=="873.c")]==1)) %in% nonsuccessful_tips))/length(nonsuccessful_tips)

pvalues_NA<-data.frame(Locus=as.numeric(gsub("\\..*","",names(out_NA$terminal$p.vals))),Significance.score=abs(out_NA$terminal$corr.dat))
pdf("Manhatten_NA.pdf",width = 6,height=5)
ggplot(pvalues_NA,aes(x=Locus,y=Significance.score))+geom_point() +
  theme_bw() + geom_hline(yintercept=as.numeric(out_NA$terminal$sig.thresh),col="red") + ylim(c(0,0.8))  +
  annotate("text",x = 240, y = 0.78, label ="240",size = 4) +
  annotate("text",x = 873, y = 0.78, label ="873",size = 4) +
  annotate("text",x = 1018, y = 0.78, label ="1018",size = 4)
dev.off()

pvalues_HA<-data.frame(Locus=as.numeric(gsub("\\..*","",names(out_HA$terminal$p.vals))),Significance.score=abs(out_HA$terminal$corr.dat))
pdf("Manhatten_hA.pdf",width = 6,height=5)
ggplot(pvalues_HA,aes(x=Locus,y=Significance.score))+geom_point() +
  theme_bw() + geom_hline(yintercept=as.numeric(out_HA$terminal$sig.thresh),col="red") + ylim(c(0,0.8))
  
dev.off()

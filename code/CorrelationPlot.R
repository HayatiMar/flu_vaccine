require(e1071)
require(corrplot)
require(ape)
require(phytools)
source("code/Tree_Statistics.R")
source("code/general_functions.R")
source("code/msvmRFE.R")

# Load data
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

#divide the data to test and train -> test on a test clade and train on other clades
test_subtree = extract.clade(tree,76876)
test_names = test_subtree$tip.label
test_ind = match(test_names,data_tip)
testi = test_ind[which(!is.na(test_ind))]
ind = which(data$outcome==1)
data$outcome[ind] = 0
data$outcome[-ind] = 1
data = scaleData(data)

train = data[-testi,]
train = train[sample(seq_len(nrow(train)), size = nrow(train)),]
train<-cbind(train$outcome,train[,1:79])

feature_selection<-svmRFE(train, k=10)
train<-train[,2:80]
write.csv(train,"trained_SVM.csv",row.names = F)
res <- cor(train)
res<-res[order(feature_selection),order(feature_selection)]
pdf("CorrelationPlot.pdf",width = 15, height = 15)
corrplot(res,type='lower',tl.pos='lt')
dev.off()

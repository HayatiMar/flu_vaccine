library(reshape2)
#library("wesanderson")
library("RColorBrewer")
library("ggplot2")
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
library(grid)


tab=read.csv("~/General_data/vaccine_subtree.csv")
colnames(tab)=c("Year","Min","Max","Median","WHO")
tab_WHO = tab[,c(1,5)]
df_p=as.data.frame(tab)
df_p.m=melt(df_p,id.vars='Year')
colnames(df_p.m)=c("Year","Approach","value")

p<-ggplot(df_p.m, aes(x=Year, y=value)) +geom_bar(aes(fill = Approach),
                                                  width = 0.2, position = position_dodge(width=0.5), stat="identity") +xlab("All Trees")+ylab("AUC")+
  ggtitle("Epitope distances between vaccine candidate sequences and \n those sequences that circulated in the following season") +
  xlab("Year") + ylab("Genetic distance")+theme_minimal()
p+theme(plot.title = element_text(size = 12, face = "bold",hjust = 0.5,lineheight = 0.9),
        legend.position="bottom", legend.title = element_blank())
#=========================================================================================
#violin plot
data = read.csv("~/General_data/var_data.csv")
data = data[,2:ncol(data)]

data$name = as.factor(data$name)
data$Model = as.factor(data$Model)
p <- ggplot()  +
  geom_violin(data = data, aes(x=name, y=value, fill=Model),width = 0.6, trim = TRUE,scale = "width",position = position_dodge(1.0))+
  geom_boxplot(data = data, aes(x=name, y=value, fill=Model), width=0.10, color="black", alpha=0.2,size = 0.3, position=position_dodge(1.0)) +
  scale_fill_viridis_d(alpha=0.7) +
  theme_minimal()+
  xlab("")+
  theme(axis.text.x = element_text(color="black", size=10, angle=0),
        plot.title = element_text(color = "black", size = 14, face = "bold",hjust = 0.5))+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14))+
  ylab("Genetic distance (amino acids)")+
  xlab("Year")
p <- p + geom_point(data = data.frame(x = factor(tab_WHO[,1]), y = tab_WHO[,2]),
                    aes(x=x, y=y), shape="square",
                    color = "#39568CFF", size = 4, alpha =0.9) +
  labs(tag = "WHO") +
  theme(plot.tag.position = c(0.947, 0.41), plot.tag = element_text(size = 9))
p
g <- rectGrob(x = 0.895, y = 0.41, width = 0.02, height = 0.03, gp = gpar(col = "#39568CFF", lty = 1, fill = "#39568CFF"))
grid.draw(g)




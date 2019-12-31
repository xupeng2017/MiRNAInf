GeneDESeq2<-function(expressioncancer,expressionnormal){
#
#miRNA/mRNA差异表达分析
print("差异表达分析DESeq2")

 
  
  
library(sqldf)
library(DESeq2)
library(tcltk)
colnames(expressioncancer)[1]<-"GeneName"
colnames(expressionnormal)[1]<-"GeneName"
expressioncancer<-sqldf("SELECT * from expressioncancer where GeneName in (select GeneName from expressionnormal) order by GeneName")
expressionnormal<-sqldf("SELECT * from expressionnormal where GeneName in (select GeneName from expressioncancer) order by GeneName")

lengthsample<-length(expressioncancer)
expre<-as.matrix(data.frame(expressioncancer[,2:lengthsample],expressionnormal[,2:lengthsample]))
expre<-round(expre)
row.names(expre)=expressioncancer[,1]

group<-matrix(seq(1,2*(lengthsample-1)))
for (i in 1:(lengthsample-1)) {
  group[i,1]<-"cancer"
  group[i+lengthsample-1,1]<-"normal"
}
group<-factor(group[,1])

colData <- data.frame(row.names=colnames(expre), group_list=group)
dds <- DESeqDataSetFromMatrix(countData = expre,
                              colData = colData,
                              design = ~ group_list)
dds2 <- DESeq(dds)
res <-  results(dds2, contrast=c("group_list","cancer","normal"))
resOrdered <- res[order(res$padj),]
resOrdered=as.data.frame(resOrdered)
resOrdered<-na.omit(resOrdered)
DEGenes<-resOrdered[resOrdered$pvalue<=0.05,]
DEGenes<-data.frame(row.names(DEGenes),DEGenes) #返回参数2
colnames(DEGenes)[1]<-"GeneName"
return(DEGenes)
}

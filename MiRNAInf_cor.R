#compute the Pearson correlation coefficient between the expression changes of a miRNA and that of its target gene

MiRNAInf_cor<-function(miRNAcancer,miRNAnormal,mRNAcancer,mRNAnormal,miRNAtargets){
  
  library(sqldf)
  library(DESeq2)
  library(tcltk)
  
  colnames(miRNAtargets)[]<-c("miRNA","Gene")
  colnames(mRNAcancer)[1]<-c("Gene")
  colnames(mRNAnormal)[1]<-c("Gene")
  
  miRNAs<-sqldf("SELECT distinct miRNA from miRNAtargets")
  mRNAs<-sqldf("SELECT distinct Gene from miRNAtargets")
  
  #去除靶基因数据表中未表达的mRNAs相关靶基因信息
  miRNAtargets1<-sqldf("SELECT * from miRNAtargets where Gene in (select Gene from mRNAnormal)")
  miRNAtargets1<-sqldf("SELECT * from miRNAtargets1 where Gene in (select Gene from mRNAcancer)")
  miRNAtargetsNUM<-sqldf("SELECT miRNA as miRNA, count(miRNA) as quantity FROM miRNAtargets1 GROUP BY miRNA")
  
  #读入miRNA前体与成熟体对照表
  preMatureMiRNA<-read.table('--/pre-mature-miRNA.txt',header = TRUE)
  colnames(preMatureMiRNA)[]<-c("premiRNA","maturemiRNA")
  
  #对miRNA表达谱进行处理
  colnames(miRNAcancer)[1]<-c("miRNA")
  miRNAcancer1<-sqldf("SELECT * from miRNAcancer,preMatureMiRNA where premiRNA=miRNAcancer.miRNA")
  miRNAcancer2<-sqldf("SELECT * from miRNAcancer1 where maturemiRNA in (select miRNA from miRNAs)")
  miRNAcancer3<-miRNAcancer2[,c(length(miRNAcancer2),2:(length(miRNAcancer2)-2))]
  miRNAcancer4<-aggregate(miRNAcancer3[,-1],data.frame(miRNAcancer3$maturemiRNA),mean)
  colnames(miRNAcancer4)[1]<-"miRNA"
  rm(miRNAcancer1,miRNAcancer2,miRNAcancer3)
  
  colnames(miRNAnormal)[1]<-c("miRNA")
  miRNAnormal1<-sqldf("SELECT * from miRNAnormal,preMatureMiRNA where premiRNA=miRNAnormal.miRNA")
  miRNAnormal2<-sqldf("SELECT * from miRNAnormal1 where maturemiRNA in (select miRNA from miRNAs)")
  miRNAnormal3<-miRNAnormal2[,c(length(miRNAnormal2),2:(length(miRNAnormal2)-2))]
  miRNAnormal4<-aggregate(miRNAnormal3[,-1],data.frame(miRNAnormal3$maturemiRNA),mean)
  colnames(miRNAnormal4)[1]<-"miRNA"
  rm(miRNAnormal1,miRNAnormal2,miRNAnormal3)
  
  #对miRNA,mRNA表达谱进行标准化（RPM）
  source('G:/xupeng/R_work/normalization.R')
  miRNAcancer4<-normalization(miRNAcancer4)
  miRNAnormal4<-normalization(miRNAnormal4)
  
  mRNAcancer<-normalization(mRNAcancer)
  mRNAnormal<-normalization(mRNAnormal)
  
  
  #miRNA表达矩阵相减
  #检验第一列名称是否完全一致
  if(all(miRNAcancer4[,1] == miRNAnormal4[,1]) ) {
    subtractionMiRNA<-as.matrix(miRNAcancer4[,-1]-miRNAnormal4[,-1])
  }else{
    miRNAcancer4<-sqldf("SELECT * from miRNAcancer4 order by miRNA")
    miRNAnormal4<-sqldf("SELECT * from miRNAnormal4 order by miRNA")
    subtractionMiRNA<-as.matrix(miRNAcancer4[,-1]-miRNAnormal4[,-1])
  }
  submiRNA<-data.frame(miRNAcancer4[,1],subtractionMiRNA)
  colnames(submiRNA)[1]<-"miRNA"
  
  #处理mRNA水平数据
  if(all(mRNAcancer[,1] == mRNAnormal[,1])) {
    subtractionMRNA<-mRNAcancer[,-1]-mRNAnormal[,-1]
  }else{
    mRNAcancer<-sqldf("SELECT * from mRNAcancer order by Gene")
    mRNAnormal<-sqldf("SELECT * from mRNAnormal order by Gene")
    subtractionMRNA<-mRNAcancer[,-1]-mRNAnormal[,-1]
    
  } 
  submRNA<-data.frame(mRNAnormal[,1],subtractionMRNA)
  colnames(submRNA)[1]<-"Gene"
  rm(subtractionMRNA)
  
  #将靶基因数据库与submRNA表格联合:miRsubmRNA
  miRsubmRNA<-sqldf("SELECT * from miRNAtargets,submRNA where submRNA.Gene=miRNAtargets.Gene")
  miRsubmRNA<-miRsubmRNA[,c(-3)]
  #将miRsubmRNA表中所有值取绝对值
  length1<-length(miRsubmRNA)
  miRsubmRNA1<-lapply(miRsubmRNA[,c(3:length1)],abs)
  #
  miRsubmRNA2<-data.frame(miRsubmRNA[,c(1,2)],miRsubmRNA1)
  rm(miRsubmRNA1)
  mRNAsum<-aggregate(miRsubmRNA2[,c(3:length1)],data.frame(miRsubmRNA2$miRNA),sum)
  colnames(mRNAsum)[1]<-"miRNA"
  
  #筛选同时出现在submiRNA以及mRNAsum表中的miRNA信息
  submiRNA1<-sqldf("SELECT * from submiRNA where submiRNA.miRNA in (select miRNA from mRNAsum) order by miRNA")
  mRNAsum1<-sqldf("SELECT * from mRNAsum where mRNAsum.miRNA in (select miRNA from submiRNA) order by miRNA")
  #submiRNA1/mRNAsum1
  divisionmiRNA<-as.matrix(submiRNA1[,-1]/mRNAsum1[,-1])
  divisionmiRNA<-data.frame(submiRNA1[,1],divisionmiRNA)
  colnames(divisionmiRNA)[1]<-"miRNA"
  
  #divisionmiRNA*miRsubmRNA2
  miRmRNAcol<-miRsubmRNA2[,c(1,2)]
  divisionmiRNA1<-sqldf("SELECT * from miRmRNAcol,divisionmiRNA where miRmRNAcol.miRNA=divisionmiRNA.miRNA order by miRNA")
  divisionmiRNA1<-divisionmiRNA1[,c(-3)]
  miRsubmRNA3<-sqldf("SELECT * from miRsubmRNA2 where miRNA in (select miRNA from submiRNA) order by miRNA")
  
  Repression<-as.matrix(divisionmiRNA1[,c(-1,-2)]*miRsubmRNA3[,c(-1,-2)])
  Repression<-data.frame(miRsubmRNA3[,c(1,2)],Repression)
  #将每个miRNA对靶基因的影响进行累加
  
  Repression1<-aggregate(Repression[,c(3:length1)],data.frame(Repression$Gene),sum)
  colnames(Repression1)[1]<-"Gene"
  Repression1<-sqldf("SELECT * from Repression1 order by Gene")
  rtRepression1<-Repression1
  
  submRNA<-sqldf("SELECT * from submRNA where submRNA.Gene in (select Gene from Repression1) order by Gene")
  
  if(all(mRNAcancer[,1] == mRNAnormal[,1])){
    cortest<-matrix(0,nrow = nrow(submRNA),ncol = 2)
    lengthsample<-length(Repression1)
    for (i in 1:nrow(submRNA)) {
      temp<-cor.test(as.numeric(Repression1[i,c(2:lengthsample)]),as.numeric(submRNA[i,c(2:lengthsample)]),alternative = "two.sided",method = "pearson",conf.level = 0.95)
      cortest[i,1]<-temp$estimate
      cortest[i,2]<-temp$p.value
    
    }
  }else{
    print("Flag1:两个表中首列基因必须一一对应")
    
  }
  
  resultcor<-data.frame(Repression1$Gene,cortest)
  colnames(resultcor)[]<-c("Gene","pearsonCor","corpvalue")
  resultcor<-na.omit(resultcor)
  resultcor<-resultcor[abs(resultcor$pearsonCor)>=0.2,]
  resultcor<-resultcor[resultcor$corpvalue<=0.05,] #返回参数1
  
  #mRNA差异表达分析
  print("差异表达分析")
  expre<-as.matrix(data.frame(mRNAcancer[,2:lengthsample],mRNAnormal[,2:lengthsample]))
  expre<-round(expre)
  row.names(expre)=mRNAcancer[,1]
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
  DEGenes<-resOrdered[resOrdered$pvalue<0.05,]
  DEGenes<-data.frame(row.names(DEGenes),DEGenes) #返回参数2
  colnames(DEGenes)[1]<-"GeneSymbol"
  
  
  #筛选miRNA-mRNA调控矩阵
  print("筛选miRNA-mRNA调控矩阵")
  #返回参数miRNAmRNAmatrix
  
  miRNAmRNAmatrix<-sqldf("SELECT * from Repression where Repression.Gene in (select GeneSymbol from DEGenes)")
  miRNAmRNAmatrix<-sqldf("SELECT * from miRNAmRNAmatrix where miRNAmRNAmatrix.Gene in (select Gene from resultcor)")
  colnames(miRNAmRNAmatrix)[c(1,2)]<-c("miRNA","Gene")
  
  #计算每个miRNA的贡献
  #返回参数miRNAcontribution
  print("计算每个miRNA的贡献")
  miRNAcontribution<-miRNAmRNAmatrix
  miRNAcontribution<-data.frame(miRNAcontribution[,c(1,2)],abs(miRNAcontribution[,c(3:length(miRNAcontribution))]))
  colnames(miRNAcontribution)[c(1,2)]<-c("miRNA","Gene")
  
  print("标记1")
  temp1<-aggregate(miRNAcontribution[,c(3:length(miRNAcontribution))],data.frame(miRNAcontribution[,2]),sum)
  colnames(temp1)[1]<-"Gene"
  
  miRNAcontributioncol<-miRNAcontribution[,c(1,2)]
  temp2<-sqldf("SELECT * from miRNAcontributioncol,temp1 where miRNAcontributioncol.Gene=temp1.Gene order by Gene")
  miRNAcontribution<-sqldf("SELECT * from miRNAcontribution order by Gene")
  
  temp2<-temp2[,c(-3)]
  #清理temp1
  rm(temp1,miRNAcontributioncol)
  
  if(all(miRNAcontribution[,2] == temp2[,2])){
    temp3<-as.matrix(miRNAcontribution[,c(-1,-2)]/temp2[,c(-1,-2)])
    miRNAcontribution<-data.frame(miRNAcontribution[,c(1,2)],temp3)
  }else{
    print("Flag2:两个表中第二列列基因必须一一对应")
    
  }
  
  
  #对miRNA表达谱进行差异表达分析
  source('G:/xupeng/R_work/GeneDESeq2.R')
  DEmiRNA<-GeneDESeq2(miRNAcancer4,miRNAnormal4)
  
  #计算miRNA对mRNA的平均影响分数Influence
  miRNAInfluence<-miRNAcontribution
  miRNAInfluence$mean<-apply(miRNAInfluence[,c(3:length(miRNAInfluence))],1,function(x){mean(x,na.rm = TRUE)})
  miRNAInfluence$NumNULL<-apply(miRNAInfluence[,c(3:(length(miRNAInfluence)-1))],1,function(x){sum(is.na(x))})
  miRNAInfluence<-miRNAInfluence[miRNAInfluence$NumNULL<((length(miRNAInfluence)-4)*0.3),c(1,2,length(miRNAInfluence)-1)]
  
  
  #将Influence表与DEmiRNA表整合
  miRNAInfluence<-sqldf("SELECT * from DEmiRNA,miRNAInfluence where DEmiRNA.GeneName=miRNAInfluence.miRNA")
  miRNAInfluence<-data.frame(miRNAInfluence$miRNA,miRNAInfluence$padj,miRNAInfluence$log2FoldChange,miRNAInfluence$mean,miRNAInfluence$Gene)
  colnames(miRNAInfluence)<-c("miRNA","miRpadj","miRlog2FC","Influence","Gene")
  
  #将Influence表与DEGenes表整合
  miRNAInfluence<-sqldf("SELECT * from miRNAInfluence,DEGenes where DEGenes.GeneSymbol=miRNAInfluence.Gene")
  miRNAInfluence<-miRNAInfluence[,c(1:5,8,12)]
  colnames(miRNAInfluence)[c(6,7)]<-c("Genelog2FC","Genepadj")
  
  #将Influence表与resultcor表整合
  miRNAInfluence<-sqldf("SELECT * from miRNAInfluence,resultcor where resultcor.Gene=miRNAInfluence.Gene")
  miRNAInfluence<-miRNAInfluence[,c(-8)]
  
  
  #读入PPI度信息表,将miRNAInfluence与PPIDegree整合
  PPIDegree<-read.csv("G:/xupeng/PPIDegree.csv",header=TRUE)
  miRNAInfluence<-sqldf("SELECT * from miRNAInfluence,PPIDegree where PPIDegree.Genesymbol=miRNAInfluence.Gene")
  miRNAInfluence<-miRNAInfluence[,c(-10)]
  
  
  
  #resultlist<-list(rtRepression1,PcormiRNAInfluence,NcormiRNAInfluence,miRNAcontribution,miRNAmRNAmatrix,submRNA,DEGenes,resultcor)
  resultlist<-list(miRNAInfluence,resultcor)
  
}
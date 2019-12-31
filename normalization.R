#标准化miRNA表达谱数据，使得表达总量为1000000；
normalization<-function(miRNA){
  for(i in 2:length(miRNA)){
    total<-sum(miRNA[,i])
    miRNA[,i]<-miRNA[,i]*10^6/total
  }
  
  return(miRNA)
  
  
  
  
  
}
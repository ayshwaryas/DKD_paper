library(tools)
library(randomForest)
library(Seurat)
library(reshape2)
library(dplyr)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(lme4)

dkd=readRDS("DKD_old_annotated_harmony_update.rds")
Idents(dkd)="celltype"
getpoissonresult <- function(celltype,cortex.integrated_sub){
  results=NULL
  gec=subset(cortex.integrated_sub,ident=celltype)
  rawdata=gec@assays$originalexp@counts[,names(gec@active.ident)]    
  rs=rawdata
  rs[rs>0]=1
  rs=as.matrix(rs)
  rsum=rowSums(rs)/(dim(rawdata)[2])
  rsum=rsum*100
  I=which(rsum>10)
  df=rawdata[I,]
  meta=data.frame(cellnames=rownames(gec@meta.data),gec@meta.data)  
  for (gene in c(1:dim(df)[1])){
    print(gene)
    data=data.frame(gene=df[gene,],meta)
    data$condition=factor(data$Condition,levels=c("non-DKD","DKD"))
    mc=mean(data$nUMI)
    data=data%>%mutate(scalenUMI=nUMI/mc)
    m <- glmer(f1, data =data, family = poisson())
    coefs <- data.frame(coef(summary(m)))
    tmp=NULL
    for(i in c(1:2)){
      if (i==1){
        x=coefs[i,]
        names(x)=paste0(row.names(coefs)[i],"_",names(x))
        tmp=data.frame(celltype=celltype, gene=rownames(df)[gene],x) #
      } else{
        x=coefs[i,]
        names(x)=paste0(row.names(coefs)[i],"_",names(x))
        tmp=data.frame(tmp,x)}
    }
    #datum=rbind(datum,data)
    results=rbind(results,tmp)
  }
  
  
  rownames(results)=results$gene
  colnames(results)=gsub("X.","",colnames(results))
  names(results)[10]="DKD_pval"
  
  return(results)
}
f1= 'gene~condition+offset(log(scalenUMI))+(1|orig.ident)'



library(tools)
library(Seurat)
library(reshape2)
library(dplyr)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(lme4)

btbr10wk_mac=readRDS("~/btbr_wk10_macs_update.rds")
btbr10wk=readRDS("~/btbr_wk10_number_annotated_analysis_update.rds")


x_lam=btbr10wk_mac@meta.data%>%group_by(orig.ident,celltype,Condition)%>%summarize(ncells=n())
tcells=btbr10wk@meta.data%>%group_by(orig.ident)%>%summarize(tcells=n())

data_btbr10wk=x_lam%>%left_join(tcells,by=c("orig.ident"))

celltypes=unique(x_lam$celltype)

f1="ncells~offset(log2(tcells))+Condition"
data_btbr10wk$Condition=factor(data_btbr10wk$Condition,levels=c("WT","OB"))
results=NULL
for (ct in celltypes[-7]){
  lam=data_btbr10wk%>%subset(celltype==ct)
  m <- glm(f1, data =lam, family = poisson())
  coefs <- data.frame(coef(summary(m)))
  tmp=NULL
  for(i in c(1:2)){
    if (i==1){
      x=coefs[i,]
      names(x)=paste0(row.names(coefs)[i],"_",names(x))
      tmp=data.frame(celltype=ct,x) #
    } else{
      x=coefs[i,]
      names(x)=paste0(row.names(coefs)[i],"_",names(x))
      tmp=data.frame(tmp,x)}
  }
  results=rbind(results,tmp)
  
}

names(results)[9]="btbr10wk_p"
names(results)[6]="btbr10wk_e"
names(results)[7]="btbr10wk_se"

padjp=p.adjust(results$btbr10wk_p,"BH")           
results=results%>%mutate(padj=padjp)

results=results%>%mutate(Significant=ifelse(padj<0.05,"Significant",""))

#### btbr10wk ####
g1<- ggplot(results, aes(x=celltype, y=btbr10wk_e, color=Significant, alpha=Significant)) + 
  geom_bar(stat="identity", width=0.4,position=position_dodge()) +geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin=btbr10wk_e-btbr10wk_se, ymax=btbr10wk_e+btbr10wk_se), width=.2,position=position_dodge(.9)) +labs(title="#cells ~ Condition", x=" ", y = "Effect-size (Condition)")+theme_classic() +scale_alpha_manual(values=c(0.75,1))+scale_color_manual(values=c("gray","black"))+ guides(fill=FALSE)+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),strip.text=element_text(size=10),axis.title.x=element_text(size=10))+coord_flip()



pdf("glm_btbr10wk.pdf",w=6,h=3)
print(g1)
dev.off()


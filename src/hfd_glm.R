
library(tools)
library(randomForest)
library(Seurat)
library(reshape2)
library(dplyr)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())
library(lme4)

hfd_mac=readRDS("~/hfd_macs_update.rds")
hfd=readRDS("~/hfd_number_annotated_analysis_update.rds")


x_lam=hfd_mac@meta.data%>%group_by(orig.ident,Age,celltype,Condition)%>%summarize(ncells=n())
tcells=hfd@meta.data%>%group_by(orig.ident,Age)%>%summarize(tcells=n())

data_hfd=x_lam%>%left_join(tcells,by=c("Age","orig.ident"))

celltypes=unique(x_lam$celltype)

f1="ncells~offset(log2(tcells))+Age+Condition"
data_hfd$Condition=factor(data_hfd$Condition,levels=c("WT","HFD"))

## removed groups with <10 cells in any condition. res mac 9, res mac 8, cycling immune


results=NULL
for (ct in celltypes[-c(1,10,12)]){
  lam=data_hfd%>%subset(celltype==ct)
  m <- glm(f1, data =lam, family = poisson())
  coefs <- data.frame(coef(summary(m)))
  tmp=NULL
  for(i in c(1:3)){
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

names(results)[9]="age_p"
names(results)[13]="hfd_p"
names(results)[6]="age_e"
names(results)[7]="age_se"
names(results)[10]="hfd_e"
names(results)[11]="hfd_se"

padjp=p.adjust(results$hfd_p,"BH")           
results=results%>%mutate(padj=padjp)

padjp=p.adjust(results$age_p,"BH")           
results=results%>%mutate(age_padj=padjp)

results=results%>%mutate(Significant=ifelse(padj<0.05,"Significant",""))

results=results%>%mutate(Significanta=ifelse(age_padj<0.05,"Significant",""))

#### HFD ####
g1<- ggplot(results, aes(x=celltype, y=hfd_e, color=Significant, alpha=Significant)) + 
  geom_bar(stat="identity", width=0.4,position=position_dodge()) +geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin=hfd_e-hfd_se, ymax=hfd_e+hfd_se), width=.2,position=position_dodge(.9)) +labs(title="#cells ~ Age + Condition", x=" ", y = "Effect-size (Condition)")+theme_classic() +scale_alpha_manual(values=c(0.75,1))+scale_color_manual(values=c("gray","black"))+ guides(fill=FALSE)+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),strip.text=element_text(size=10),axis.title.x=element_text(size=10))+coord_flip()


g2<- ggplot(results, aes(x=celltype, y=age_e, color=Significanta, alpha=Significanta)) + 
  geom_bar(stat="identity", width=0.4,position=position_dodge()) +geom_hline(yintercept = 0)+
  geom_errorbar(aes(ymin=age_e-age_se, ymax=age_e+age_se), width=.2,position=position_dodge(.9)) +labs(title="#cells ~ Age + Condition", x=" ", y = "Effect-size (Age)")+theme_classic() +scale_alpha_manual(values=c(0.75,1))+scale_color_manual(values=c("gray","black"))+ guides(fill=FALSE)+theme(axis.text.x=element_text(size=10),axis.text.y=element_text(size=10),strip.text=element_text(size=10),axis.title.x=element_text(size=10))+coord_flip()

g=plot_grid(g1,g2,ncol=2)

pdf(paste0("glmhfd.pdf"),w=12,h=3)
print(g)
dev.off()



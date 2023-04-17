library(tools)
library(randomForest)
library(Seurat)
library(reshape2)
library(dplyr)
library(tidyverse)
library(cowplot)
theme_set(theme_cowplot())

b10_mac=readRDS("btbr_wk10_macs_update.rds")
b5_mac=readRDS("/btbr_wk5_macs_update.rds")

b10=FindVariableFeatures(b10_mac, selection.method = "vst",nfeatures = 3000)
b4=FindVariableFeatures(b5_mac, selection.method = "vst",nfeatures = 3000)

genes.use=intersect(b4@assays$RNA@var.features,b10@assays$RNA@var.features) #1489
training.set = c()
test.set=c()
training.label = c()
test.label=c()

for (i in (levels(b10_mac@active.ident))){
  cells.in.clust = WhichCells(b10_mac,idents =i);
  n = min(1000, round(length(cells.in.clust)*0.7))
  train.temp = cells.in.clust[sample(length(cells.in.clust))][1:n]
  test.temp = setdiff(cells.in.clust, train.temp)
  training.set = c(training.set,train.temp); test.set=c(test.set,test.temp)
  training.label = c(training.label, rep(i,length(train.temp))); test.label = c(test.label, rep(i, length(test.temp)));
}

predictor_Data = as.matrix(b10_mac@assays$RNA@data[genes.use,])
tmp = as.vector(table(training.label))
sampsizes = rep(min(tmp),length(tmp))
predictor_Data = t(scale(t(predictor_Data), center=TRUE, scale=TRUE))
rf_b10_mac = randomForest(x=t(predictor_Data[,training.set]), y=factor(training.label), importance = TRUE, ntree = 5001, proximity=TRUE, sampsize=sampsizes, keep.inbag=TRUE, replace=FALSE)


save(rf_b10_mac,file="updated_rf_b10_mac.RData")

test.predict = predict(rf_b10_mac,t(predictor_Data[,test.set]))
Conf_test = table(test.label,test.predict)
save(Conf_test,genes.use,file="rf_b10_mac_trainset.RObj")

###### predict on b4 ######

hca.rf = b5_mac@assays$RNA@data[genes.use, ]
hca.rf = t(scale(t(as.matrix(hca.rf)), center=TRUE, scale=TRUE))
hca.rf[is.na(hca.rf)] = 0
hca.predict = predict(rf_b10_mac,t(hca.rf))

hca.ident = factor(b5_mac@active.ident)#factor(df1$Cluster) #factor(pooled_kidney@active.ident)
Conf_test = table(hca.ident, hca.predict[names(hca.ident)])
Conf_test_M=melt(Conf_test,id=rownames(Conf_test))
Conf_test_M_sum=data.frame(hca.ident=rownames(Conf_test),sums=rowSums(Conf_test))
Conf_test_M=Conf_test_M%>%left_join(Conf_test_M_sum,by="hca.ident")
Conf_test_M=Conf_test_M%>%mutate(prop=(value/sums))

names(Conf_test_M)[2]="ident"


Conf_test_M$hca.ident=factor(Conf_test_M$hca.ident,levels=c("Infiltrating Macrophages","Resident Macrophage-1","Resident Macrophage-2","Resident Macrophage-7","Resident Macrophage-8/Doublet","Proliferating macrophages"))
Conf_test_M$ident=factor(Conf_test_M$ident,levels=c("Infiltrating Macrophage-1","Infiltrating Macrophage-2","Resident Macrophage-1","Resident Macrophage-2","Resident Macrophage-3","Resident Macrophage-5","Resident Macrophage-6","Proliferating macrophages"))


### plotting results


pdf("rf_b10_mac_b5.pdf",8,6,useDingbats = FALSE)
Conf_test_M%>%ggplot(aes(x=ident,y=hca.ident))+xlab("Predicted cell-types")+ylab("Annotation")+geom_point(aes(size=(prop*100+1),color=(prop*100+1)))+theme(axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5))+ scale_color_gradient(low="white", high="red")+guides(size=guide_legend(title="Proportion of cells"),color=guide_legend(title=""))
dev.off()


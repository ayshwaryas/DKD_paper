
library(pheatmap)
library(tidyr)

generate_gene_df <- function(seurat_obj,genes,df1){
  I=which(rownames(seurat_obj@assays$RNA@data)%in%genes$Genes)
  df=seurat_obj@assays$RNA@data[I,]
  df_m=melt(as.matrix(df))
  df1=data.frame(Var2=colnames(seurat_obj@assays$RNA@data),df1)
  head(df1)
  df_m=df_m%>%left_join(df1,by="Var2")
  names(df_m)[1]="genes"
  df_m$genes <- factor(df_m$genes, levels = rev(genes$Genes))
  return(df_m)
}

tf3=data.frame(Genes=c("C1qa","C1qb","C1qc","Cd74","Fcgr1","Cd81","ApoE","Fyb","Adgre1","Ms4a7","Lgmn","Hexb","Cd72","Cxcl16","Ccl12","Slamf9","Ccr2","Chil3","Fabp4","Ear2","Crip1","S100a11","S100a4","S100a6","Napsa","Gpx1","Vim","Smpdl3a","Ifitm6","Plac8","Adgre4","Nupr1"))

b10_mac=readRDS("btbr_wk10_macs_update.rds")
b5_mac=readRDS("btbr_wk5_macs_update.rds")


df_hfd=hfd_mac@meta.data
df_b10=b10_mac@meta.data
df_b5=b5_mac@meta.data



#### TF3
hfd_genes=generate_gene_df(hfd_mac,tf3,df_hfd)
b10_genes=generate_gene_df(b10_mac,tf3,df_b10)
b4_genes=generate_gene_df(b5_mac,tf3,df_b5)


data=data.frame(group="HFD",hfd_genes%>%select(Var2,genes,value,celltype))
data=rbind(data,data.frame(group="BTBR10wk",b10_genes%>%select(Var2,genes,value,celltype)))
data=rbind(data,data.frame(group="BTBR5wk",b4_genes%>%select(Var2,genes,value,celltype)))

data$group=factor(data$group,levels=c("BTBR5wk","BTBR10wk","HFD"))

data$celltype=gsub(" ","",data$celltype)
data$celltype=gsub("-","",data$celltype)
data$celltype=gsub("_","",data$celltype)

df_m=data%>%group_by(genes,celltype,group)%>%summarise(avgexp=mean(value))



df_m$celltype=factor(df_m$celltype,levels=unique(df_m$celltype))
df_m$group=factor(df_m$group,levels=c("BTBR5wk","BTBR10wk","HFD"))
df_m=df_m%>%arrange(group,celltype)
df_mu=unite(df_m,col="celltype_group",celltype,group)
df_mu$celltype_group=factor(df_mu$celltype_group,levels=unique(df_mu$celltype_group))
df_mu_s=spread(df_mu,celltype_group,avgexp)
mat=apply(df_mu_s[,-1],1:2,as.numeric)
rownames(mat)=df_mu_s$genes
mat_col <- data.frame(group =sapply(colnames(df_mu_s)[-1],function(x) return(strsplit(x,"_")[[1]][2])))
mat_row <- data.frame(Gene = df_mu_s$genes)
mat_comp=data.frame(celltype =sapply(colnames(df_mu_s)[-1],function(x) return(strsplit(x,"_")[[1]][1])))

ann.col <- list(group = c("BTBR5wk"="yellow","BTBR10wk"="purple","HFD"="gray"))
pheatmap(mat,filename="dkd_mm_immmunegenes_heatmap_residentinf.pdf",color = colorRampPalette(c("blue", "white", "red"))(13),scale='row',cluster_rows = FALSE, cluster_cols = FALSE, annotation_colors = ann.col,labels_col = as.character(mat_comp$celltype), width = 8, height = 8, fontsize = 10, angle_col=90, annotation_col = mat_col)






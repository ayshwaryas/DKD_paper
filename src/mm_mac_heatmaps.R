library(tidyr)
library(pheatmap)

btbr_5wk_markers=FindAllMarkers(b5_mac)
btbr_10wk_markers=FindAllMarkers(b10_mac)
hfd_markers=FindAllMarkers(hfd_mac)

write_tsv(hfd_markers,"hfd_mac_markers.tsv")
write_tsv(btbr_5wk_markers,"btbr_5wk_mac_markers.tsv")
write_tsv(btbr_10wk_markers,"/btbr_10wk_mac_markers.tsv")

y=btbr_5wk_markers%>%subset(cluster=="Resident Macrophage-1")%>%arrange(desc(avg_log2FC))



###### Heatmap BTBR 5wk########

genes=read_tsv("btbr5wk_markers.txt")
data=generate_gene_df(b5_mac,genes,df_b5)
df_m=data%>%group_by(genes,celltype)%>%summarise(avgexp=mean(value))
df_m$celltype=factor(df_m$celltype,levels=c("Inf Macrophage","Resident Macrophage-1","Resident Macrophage-2","Resident Macrophage-7","Resident Macrophage-8/Doublet","Proliferating macrophages"))
df_m$genes=factor(df_m$genes,levels=unique(genes$Genes))
df_mu_s=spread(df_m,celltype,avgexp)
mat=apply(df_mu_s[,-1],1:2,as.numeric)
rownames(mat)=df_mu_s$genes
pheatmap(t(mat),filename="btbr5wk_immune_heatmap.pdf",color = colorRampPalette(c("blue","white", "red"))(13),scale='column',cluster_rows = FALSE, cluster_cols = FALSE, width = 12, height = 2, fontsize = 10, angle_col=45)

###### Heatmap BTBR 10wk########

genes=read_tsv("btbr10wk_markers.txt")
data=generate_gene_df(b10_mac,genes,df_b10)
df_m=data%>%group_by(genes,celltype)%>%summarise(avgexp=mean(value))
df_m$celltype=factor(df_m$celltype,levels=c("Infiltrating Macrophage-1","Infiltrating Macrophage-2","Resident Macrophage-1","Resident Macrophage-2","Resident Macrophage-3","Resident Macrophage-5","Proliferating macrophages"))
df_m$genes=factor(df_m$genes,levels=unique(genes$Genes))
df_mu_s=spread(df_m,celltype,avgexp)
mat=apply(df_mu_s[,-1],1:2,as.numeric)
rownames(mat)=df_mu_s$genes
pheatmap(t(mat),filename="btbr10wk_immune_heatmap.pdf",color = colorRampPalette(c("blue","white", "red"))(13),scale='column',cluster_rows = FALSE, cluster_cols = FALSE, width = 12, height = 2, fontsize = 10, angle_col=45)


###### Heatmap HFD########

genes=read_tsv("/hfd_markers.txt")
data=generate_gene_df(hfd_mac,genes,df_hfd)
df_m=data%>%group_by(genes,celltype)%>%summarise(avgexp=mean(value))
df_m$celltype=factor(df_m$celltype,levels=c("Infiltrating Macrophage-1","Infiltrating Macrophage-2","Infiltrating Macrophage-3","Resident Macrophage-1","Resident Macrophage-2","Resident Macrophage-3","Resident Macrophage-4","Resident Macrophage-5","Resident Macrophage-6","Resident Macrophage-8/Doublet","Resident Macrophage-9","Cycling Immune"))
df_m$genes=factor(df_m$genes,levels=unique(genes$Genes))
df_mu_s=spread(df_m,celltype,avgexp)
mat=apply(df_mu_s[,-1],1:2,as.numeric)
rownames(mat)=df_mu_s$genes
pheatmap(t(mat),filename="hfd_immune_heatmap.pdf",color = colorRampPalette(c("blue","white", "red"))(13),scale='column',cluster_rows = FALSE, cluster_cols = FALSE, width = 13, height = 3, fontsize = 10, angle_col=45)

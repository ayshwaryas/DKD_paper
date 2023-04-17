





############# VOLCANO PLOT ##############
library(ggrepel)

lam_markers=FindMarkers(hfd_mac,"Resident Macrophage-6",logfc.threshold = -10)
lam_markers=data.frame(genes=rownames(lam_markers),lam_markers)
saveRDS(lam_markers,"lam_markers.rds" )

LAM=c("Lpl","Cd5l","Trem2","Slc40a1","Fabp5","Gdf15","Lgals3","Ctsd","Ctsb","Cd9","Apoe","Vcam1","Cd38","Fcna","Lipa","Mgl2","Cd14","Il1b","Cst3","Hexb","Cd74","H2-Eb1","H2-Aa","H2-Ab1","Tmem176a","Ccr2","Lgals1","Clec1b","Clec4n","Fcna")


dsub=lam_markers
dsub=dsub%>%mutate(label=ifelse(genes%in%LAM,as.character(genes),""))
dsub=dsub%>%mutate(sig=ifelse(p_val_adj<0.05,"Adjusted p-value < 0.05","Adjusted p-value >= 0.05"))
I=which(dsub$p_val_adj==0)
dsub$p_val_adj[I]=1e-80

p = ggplot(dsub, aes(avg_log2FC, -log10(p_val_adj))) + geom_point(aes(col=sig),size=0.75) +xlab("Avg log2 Fold Change") + ylab("-log10 Adjusted p-value")
p=p+geom_text_repel(data=dplyr::filter(dsub, p_val_adj<0.05), aes(label=label),size=6)+theme(legend.title=element_blank())+scale_color_manual(values=c("dodgerblue","gray"))
pdf("LAM-volcano_wilcox.pdf",w=10,h=7,useDingbats = FALSE)
print(p)
dev.off()

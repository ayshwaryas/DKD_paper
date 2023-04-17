


b10_mac=readRDS("btbr_wk10_macs_update.rds")
b5_mac=readRDS("btbr_wk5_macs_update.rds")

df_b10=data.frame(b10_mac@reductions$phate@cell.embeddings,b10_mac@meta.data)
df_b10$celltype=factor(df_b10$celltype,levels=c("Resident Macrophage-1","Resident Macrophage-2","Resident Macrophage-3","Resident Macrophage-5","Infiltrating Macrophage-1","Infiltrating Macrophage-2","Proliferating macrophages"))
g1=df_b10%>%ggplot(aes(x=phate_1,y=phate_2,col=celltype))+geom_point(size=2)+ guides(colour = guide_legend(ncol=1,override.aes = list(size=10),title=""))+ theme(axis.title=element_text(size=16),axis.text.x = element_text(size=16),axis.text.y=element_text(size=16),legend.text=element_text(size=10),legend.title=element_blank(),strip.background = element_rect(fill="white"),strip.text.x = element_text(size = 20, colour = "black", angle = 90)) +scale_colour_manual(values=c("#F8C4B4","#FC4E2A","darkslategrey","cadetblue1","dodgerblue2","#a26d41","#665A48"))                                                                                                                                                                                                                                                                           
g2=df_b10%>%ggplot(aes(x=phate_1,y=phate_2,col=scDblFinder.class))+geom_point(size=2)+ guides(colour = guide_legend(ncol=1,override.aes = list(size=10),title=""))+ theme(axis.title=element_text(size=16),axis.text.x = element_text(size=16),axis.text.y=element_text(size=16),legend.text=element_text(size=10),legend.title=element_blank(),strip.background = element_rect(fill="white"),strip.text.x = element_text(size = 20, colour = "black", angle = 90)) +scale_colour_manual(values=c("brown","salmon"))

g=plot_grid(g1,g2,ncol=2,rel_widths = c(1.2,1))
pdf("/Volumes/broad_grekalab/ayshwarya/2023-04-02.b10_mac_phate_doublet.pdf",width=12,height=5,useDingbats=FALSE)
g
dev.off()



df_b5=data.frame(b5_mac@reductions$phate@cell.embeddings,b5_mac@meta.data)

df_b5$celltype=factor(df_b5$celltype,levels=c("Resident Macrophage-1","Resident Macrophage-2","Resident Macrophage-7","Resident Macrophage-8/Doublet","Inf Macrophage","Proliferating macrophages"))

g1=df_b5%>%ggplot(aes(x=phate_1,y=phate_2,col=celltype))+geom_point(size=2)+ guides(colour = guide_legend(ncol=1,override.aes = list(size=10),title=""))+ theme(axis.title=element_text(size=16),axis.text.x = element_text(size=16),axis.text.y=element_text(size=16),legend.text=element_text(size=10),legend.title=element_blank(),strip.background = element_rect(fill="white"),strip.text.x = element_text(size = 20, colour = "black", angle = 90)) +scale_colour_manual(values=c("#F8C4B4","#FC4E2A","salmon","#D6604D","blanchedalmond","#665A48"))
g2=df_b5%>%ggplot(aes(x=phate_1,y=phate_2,col=scDblFinder.class))+geom_point(size=2)+ guides(colour = guide_legend(ncol=1,override.aes = list(size=10),title=""))+ theme(axis.title=element_text(size=16),axis.text.x = element_text(size=16),axis.text.y=element_text(size=16),legend.text=element_text(size=10),legend.title=element_blank(),strip.background = element_rect(fill="white"),strip.text.x = element_text(size = 20, colour = "black", angle = 90)) +scale_colour_manual(values=c("brown","salmon"))

g=plot_grid(g1,g2,ncol=2,rel_widths = c(1.2,1))
pdf("b5_mac_phate_doublet.pdf",width=12,height=5,useDingbats=FALSE)
g
dev.off()


#### HFD
hfd_mac=readRDS("hfd_macs_update.rds")
df_mac=data.frame(hfd_mac@reductions$phate@cell.embeddings,hfd_mac@meta.data)
df_mac$celltype=factor(df_mac$celltype,levels=c("Resident Macrophage-1","Resident Macrophage-2","Resident Macrophage-3","Resident Macrophage-4","Resident Macrophage-5","Resident Macrophage-6","Resident Macrophage-8/Doublet","Resident Macrophage-9","Infiltrating Macrophage-1","Infiltrating Macrophage-2","Infiltrating Macrophage-3","Cycling Immune"))

df_mac$Condition=factor(df_mac$Condition,levels=c("WT","HFD"))

df_mac$Age=factor(df_mac$Age,levels=c(66,80,83,94,100))


g1=df_mac%>%ggplot(aes(x=phate_1,y=phate_2,col=celltype))+geom_point(size=2)+ guides(colour = guide_legend(ncol=1,override.aes = list(size=10),title=""))+ theme(axis.title=element_text(size=16),axis.text.x = element_text(size=16),axis.text.y=element_text(size=16),legend.text=element_text(size=10),legend.title=element_blank(),strip.background = element_rect(fill="white"),strip.text.x = element_text(size = 20, colour = "black", angle = 90)) +scale_colour_manual(values=c("#F8C4B4","#FC4E2A","darkslategrey","darkcyan","cadetblue1","#D5CEA3","#D6604D","#053061","dodgerblue2","#a26d41","#9EC9E2","#F4BFBF"))
g2=df_mac%>%ggplot(aes(x=phate_1,y=phate_2,col=scDblFinder.class))+geom_point(size=2)+ guides(colour = guide_legend(ncol=1,override.aes = list(size=10),title=""))+ theme(axis.title=element_text(size=16),axis.text.x = element_text(size=16),axis.text.y=element_text(size=16),legend.text=element_text(size=10),legend.title=element_blank(),strip.background = element_rect(fill="white"),strip.text.x = element_text(size = 20, colour = "black", angle = 90)) +scale_colour_manual(values=c("brown","salmon"))

g3=df_mac%>%ggplot(aes(x=phate_1,y=phate_2,col=Age))+geom_point(size=2)+ guides(colour = guide_legend(ncol=1,override.aes = list(size=10),title=""))+ theme(axis.title=element_text(size=16),axis.text.x = element_text(size=16),axis.text.y=element_text(size=16),legend.text=element_text(size=10),legend.title=element_blank(),strip.background = element_rect(fill="white"),strip.text.x = element_text(size = 20, colour = "black", angle = 90)) +scale_colour_manual(values=c("#9F73AB","#F1A661","#BAD1C2","#DFD3C3","#50577A"))

g4=df_mac%>%ggplot(aes(x=phate_1,y=phate_2,col=Condition))+geom_point(size=2)+ guides(colour = guide_legend(ncol=1,override.aes = list(size=10),title=""))+ theme(axis.title=element_text(size=16),axis.text.x = element_text(size=16),axis.text.y=element_text(size=16),legend.text=element_text(size=10),legend.title=element_blank(),strip.background = element_rect(fill="white"),strip.text.x = element_text(size = 20, colour = "black", angle = 90)) +scale_colour_manual(values=c("#112B3C","#FBC252"))

g=plot_grid(g1,g2,g4,ncol=3,rel_widths = c(1.3,1,1))
pdf("hfd_mac_phate.pdf",width=22,height=5,useDingbats=FALSE)
g
dev.off()


pdf("hfd_mac_phate_doublet.pdf",width=6,height=5,useDingbats=FALSE)
g2
dev.off()
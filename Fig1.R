library(Matrix)
library(Seurat)
library(tidyverse)
library(liger)
library(patchwork)
library(cowplot)
library(patchwork)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(DoubletFinder)
library(Matrix)
library(magrittr)
library(stringr)
library(tidyverse)
library(scRNAtoolVis)

load("AB_RENAME.RData")

#---------UMAP(celltype)
p7 <- DimPlot(AB, reduction = "umap", group.by = "age", pt.size=0.01,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap_celltype.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

#-----------celltype-Dotplot
jjplot <- jjDotPlot(object = AB,
                    gene = c("Rbfox3","Satb2","Slc17a7","Gad1","Gad2","Reln","Aqp4","Gfap","Slc1a3","Pdgfra","Cspg4","Mbp","Bcas1","Plp1","Csf1r","Lgmn","P2ry12","Igfbp7","Adgrl4","Vwf","Pdgfrb","Col1a1","Col1a2"),
                    gene.order = c("Rbfox3","Satb2","Slc17a7","Gad1","Gad2","Reln","Aqp4","Gfap","Slc1a3","Pdgfra","Cspg4","Mbp","Bcas1","Plp1","Csf1r","Lgmn","P2ry12","Igfbp7","Adgrl4","Vwf","Pdgfrb","Col1a1","Col1a2"),
                    id = 'celltype',
                    cluster.order = c("ExN",	"InN",	"CR",	"AST",	"OPC",	"OL",	"MG",	"Endo",	"Peri",	"VLMC"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)#+ggplot2:::coord_flip()
jjplot

ggsave(filename = "Dotplot-all.pdf", plot = jjplot, device = 'pdf', width = 25, height = 20, units = 'cm')


#-----------cellcycle
library(Seurat)

cc.genes=readxl::read_xlsx("cell_cycle.xlsx") #in supplementary material

stringr::str_to_title(cc.genes)

g1s.genes <- na.omit(cc.genes$G1_S_gene)
g2m.genes <- cc.genes$G2_M_gene

g1s.genes <- stringr::str_to_title(g1s.genes)
g2m.genes <- stringr::str_to_title(g2m.genes)

DefaultAssay(AB)="RNA"
AB <- CellCycleScoring(AB, s.features = g1s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(AB[[]])

#
RidgePlot(AB, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)

table(AB$celltype,AB$Phase)
table(AB$celltype,AB$Phase)

table(AB$celltype)#
prop.table(table(Idents(AB)))
table(Idents(AB), AB$celltype)#
Cellratio <- prop.table(table(Idents(AB), AB$celltype), margin = 2)#
Cellratio <- as.data.frame(Cellratio)
#Cellratio$Var2<-factor(Cellratio$Var2, levels=c("Normal","BCN","BCL","BCH"))
allcolour <-  c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF",
                "#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080",
                "#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD",
                "#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C",
                "#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
my36colors <- c('#d86967', '#bbbbd6', '#58539f', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#000000')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = my36colors)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),axis.text.x = element_text(size = rel(1), angle = 45, vjust = 1, hjust = 1))
ggsave(filename = "celltype_Group_Cell_cycle proportion.pdf", device = 'pdf', width = 13, height = 6, units = 'cm')

#-----------density plot
KS_plot_density(obj=AB, 
                marker= c("Top2a","Cenpf"),
                dim = "UMAP", size =0.5)

P1 <- KS_plot_density(obj=AB, 
                      marker=c("Mki67", "Cenpf","Kif20b", "Hjurp", "Kif2c","Cdk1"), 
                      dim = "UMAP", size =0.5, ncol = 3)

ggsave('density-UMAP.pdf',plot = P1, width = 18,height = 12)


#时间点DEGs
# 寻找差异表达基因

AB$celltype.age <- paste(Idents(AB), AB$age, sep = "_")
AB$celltype <- Idents(AB)
Idents(AB) <- "celltype.age"
###--
#Neuron
age.response <- FindMarkers(AB, ident.1 = c("ExN_3d"), ident.2 = ("ExN_Sham3d"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="ExN_3d_Sham3d.csv")
age.response <- FindMarkers(AB, ident.1 = c("ExN_3d"), ident.2 = ("ExN_12h"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="ExN_3d_12h.csv")
age.response <- FindMarkers(AB, ident.1 = c("ExN_3d"), ident.2 = ("ExN_7d"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="ExN_3d_7d.csv")




age.response <- FindMarkers(AB, ident.1 = c("InN_3d"), ident.2 = ("InN_Sham3d"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="InN_3d_Sham3d.csv")
age.response <- FindMarkers(AB, ident.1 = c("InN_3d"), ident.2 = ("InN_12h"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="InN_3d_12h.csv")
age.response <- FindMarkers(AB, ident.1 = c("InN_3d"), ident.2 = ("InN_7d"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="InN_3d_7d.csv")



#Glia
age.response <- FindMarkers(AB, ident.1 = c("AST_3d"), ident.2 = ("AST_Sham3d"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="AST_3d_Sham3d.csv")
age.response <- FindMarkers(AB, ident.1 = c("AST_3d"), ident.2 = ("AST_12h"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="AST_3d_12h.csv")
age.response <- FindMarkers(AB, ident.1 = c("AST_3d"), ident.2 = ("AST_7d"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="AST_3d_7d.csv")


age.response <- FindMarkers(AB, ident.1 = c("MG_3d"), ident.2 = ("MG_Sham3d"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="MG_3d_Sham3d.csv")
age.response <- FindMarkers(AB, ident.1 = c("MG_3d"), ident.2 = ("MG_12h"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="MG_3d_12h.csv")
age.response <- FindMarkers(AB, ident.1 = c("MG_7d"), ident.2 = ("MG_Sham7d"),verbose = FALSE)
head(age.response, n = 30)
write.csv(age.response,file="MG_7d_Sham7d.csv")


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

#——————————Data QC——————————#
#delete Low-quality and RBC and Doubles
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("1",	"2",	"3",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20",	"21",	"22",	"23",	"24",	"25",	"26",	"27",	"29",	"30",	"31",	"32",	"33",	"34",	"35",	"36", "37", "38", "39",	"40",	"42", "43")]

##----------------------
#####
###delete Hb+ gene
count <- AB@assays$RNA@counts
count[1:5,1:5]
count2 <- -which(rownames(count) %in% c("Hba-a2","Hba-a3","Hbb-b1","Hbb-bs","Hbg1","Hbe1","Hbe2","Hbz"))
count3 <- count[-c(2295,2298,12209,12222,19194),]

meta <- AB@meta.data

AB <- CreateSeuratObject(counts = count3, meta.data = meta , project = "HIE", min.cells = 3, min.features = 200)
###
###
myfunction1 <- function(xxxx){
  xxxx <- NormalizeData(xxxx, normalization.method = "LogNormalize", scale.factor = 10000)
  xxxx <- FindVariableFeatures(xxxx, selection.method = "vst", nfeatures = 2000)
  return(xxxx)   }

AB <- myfunction1(AB)

object.list <- SplitObject(AB, split.by = "orig.ident")

AB.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = 2000,dims = 1:30)  ##耗时久
AB <- IntegrateData(anchorset = AB.anchors, dims = 1:30)
###
DefaultAssay(AB) <- "integrated"
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot <- ElbowPlot(AB)
ggsave(filename = "ElbowPlot.pdf", plot = plot, device = 'pdf', width = 14, height = 12, units = 'cm')

save(AB,file="AB2.RData")
#
AB <- FindNeighbors(AB, dims = 1:30)
AB <- FindClusters(AB, resolution = 1.2)
AB <- RunUMAP(AB, dims = 1:30)
DimPlot(AB, reduction = "umap")
#AB <- RunTSNE(AB, dims = 1:30)
#DimPlot(AB, reduction = "tsne")
#
AB@meta.data[["age"]]<-factor(AB@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
DefaultAssay(AB) <- "RNA"
table(AB@meta.data[["seurat_clusters"]], AB@meta.data[["age"]])

save(AB,file="AB2_new.RData")

p7 <- DimPlot(AB, reduction = "umap", group.by = "age", pt.size=0.01)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2_Cluster2.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

jjplot <- jjDotPlot(object = AB,
                    gene = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1","Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Nrp2","Ntf3","Ntng1"),
                    gene.order = c("Hba-a2","Hbb-b1","Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1","Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Nrp2","Ntf3","Ntng1"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot2.pdf", plot = jjplot, device = 'pdf', width = 40, height = 40, units = 'cm')


AB100=subset(AB,downsample=100)
markers100 <- FindAllMarkers(AB100, logfc.threshold = 0.25, min.pct = 0.5, only.pos = F, test.use = "wilcox") 
write.table(markers100,file="markers100-2.txt",quote=F,sep="\t",row.names=F,col.names=T)

#——————————Data QC——————————#
#deleteLow-quality and RBC and Doubles
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20",	"21",	"22",	"23",	"24",	"26",	"27",	"29",	"30",	"31",	"32",	"33",	"34",	"35",	"36", "37", "38", "39",	"40",	"42", "43")]


###
count <- AB@assays$RNA@counts
#
meta <- AB@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","age")]
#
AB <- CreateSeuratObject(counts = count, meta.data = meta , project = "HIE", min.cells = 3, min.features = 200)
###
###
myfunction1 <- function(xxxx){
  xxxx <- NormalizeData(xxxx, normalization.method = "LogNormalize", scale.factor = 10000)#
  xxxx <- FindVariableFeatures(xxxx, selection.method = "vst", nfeatures = 2000)#
  return(xxxx)   }

AB <- myfunction1(AB)

###
object.list <- SplitObject(AB, split.by = "orig.ident")#

AB.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = 2000,dims = 1:30)
AB <- IntegrateData(anchorset = AB.anchors, dims = 1:30)
#
DefaultAssay(AB) <- "integrated"
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot <- ElbowPlot(AB)
ggsave(filename = "ElbowPlot.pdf", plot = plot, device = 'pdf', width = 14, height = 12, units = 'cm')

save(AB,file="AB3.RData")
#
AB <- FindNeighbors(AB, dims = 1:22)
AB <- FindClusters(AB, resolution = 1.2)
AB <- RunUMAP(AB, dims = 1:22)
DimPlot(AB, reduction = "umap")
AB <- RunTSNE(AB, dims = 1:22)
DimPlot(AB, reduction = "tsne")
#
AB@meta.data[["age"]]<-factor(AB@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
DefaultAssay(AB) <- "RNA"
table(AB@meta.data[["celltype"]], AB@meta.data[["age"]])

save(AB,file="AB3_new.RData")

p7 <- DimPlot(AB, reduction = "umap", group.by = "age", pt.size=0.01)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2_Cluster3.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

jjplot <- jjDotPlot(object = AB,
                    gene = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1","Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Nrp2","Ntf3","Ntng1","Nos1","Dcx","Meis2","Ndnf"),
                    gene.order = c("Hba-a2","Hbb-b1","Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1","Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Nrp2","Ntf3","Ntng1","Nos1","Dcx","Meis2","Ndnf"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot3.pdf", plot = jjplot, device = 'pdf', width = 40, height = 40, units = 'cm')

AB100=subset(AB,downsample=100)
markers100 <- FindAllMarkers(AB100, logfc.threshold = 0.25, min.pct = 0.5, only.pos = F, test.use = "wilcox")  #
write.table(markers100,file="markers100-3.txt",quote=F,sep="\t",row.names=F,col.names=T)



#——————————Renaming the cluster——————————#
#
###
current.cluster.ids <- c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20",	"21",	"22",	"23",	"24",	"25",	"26",	"27",	"28",	"29",	"30",	"31",	"32",	"33",	"34",	"35",	"36")
new.cluster.ids <- c("AST",	"ExN",	"ExN",	"OPC",	"ExN",	"ExN",	"AST",	"ExN",	"ExN",	"InN",	"ExN",	"ExN",	"MG",	"InN",	"ExN",	"InN",	"ExN",	"ExN",	"InN",	"OL",	"ExN",	"InN",	"Endo",	"ExN",	"ExN",	"AST",	"ExN",	"VLMC",	"OPC",	"ExN",	"ExN",	"VLMC",	"Peri",	"MG",	"CR",	"Endo",	"AST") 

AB@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(AB@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(AB@meta.data$celltype)

#
AB <- RenameIdents(AB, `0`="AST",	`1`="ExN",	`2`="ExN",	`3`="OPC",	`4`="ExN",	`5`="ExN",	`6`="AST",	`7`="ExN",	`8`="ExN",	`9`="InN",	`10`="ExN",	`11`="ExN",	`12`="MG",	`13`="InN",	`14`="ExN",	`15`="InN",	`16`="ExN",	`17`="ExN",	`18`="InN",	`19`="OL",	`20`="ExN",	`21`="InN",	`22`="Endo",	`23`="ExN",	`24`="ExN",	`25`="AST",	`26`="ExN",	`27`="VLMC",	`28`="OPC",	`29`="ExN",	`30`="ExN",	`31`="VLMC",	`32`="Peri",	`33`="MG",	`34`="CR",	`35`="Endo",	`36`="AST")

Idents(AB) <- factor(Idents(AB), levels = c("ExN",	"InN",	"CR",	"AST",	"OPC",	"OL",	"MG",	"Endo",	"Peri",	"VLMC"))
DimPlot(AB, reduction = "umap", label = TRUE)

table(AB@meta.data[["celltype"]], AB@meta.data[["age"]])
#####
#
save(AB,file="AB3.RData")

DefaultAssay(AB) <- "RNA"
AB100=subset(AB,downsample=100)
markers <- FindAllMarkers(AB, logfc.threshold = 0.25, min.pct = 0.5, only.pos = F, test.use = "wilcox")
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)

#---------------conjoint analysis---------------#
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
library(homologene)

#-------------MG7-------------#
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^Mt-")
###
myfunction1 <- function(xxxx){
  xxxx <- NormalizeData(xxxx, normalization.method = "LogNormalize", scale.factor = 10000)
  xxxx <- FindVariableFeatures(xxxx, selection.method = "vst", nfeatures = 2000)
  return(xxxx)   }
AB <- scRNA
AB <- myfunction1(AB)
#
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
plot2
save(AB,file="AB.RData")
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:17,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:17)
AB <- FindClusters(AB, resolution = 0.4)
AB$group=str_replace(AB$age,"Sham.*","Sham")
#
AB@meta.data[["age"]]<-factor(AB@meta.data[["age"]], levels=c("Sham","S1d","S3d","S7d"))
#
DefaultAssay(AB) <- "RNA"
AB$age <- as.character(AB@meta.data[["age"]])
save(AB,file="AB_new.RData")
#--
p7 <- DimPlot(AB, reduction = "umap", group.by = "group", pt.size=1,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=1, label = TRUE,repel = TRUE,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2_Cluster2.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

jjplot <- jjDotPlot(object = AB,
                    gene = c("Csf1r","P2ry12","Tmem119","Itgax","Spp1","Gpnmb","Piezo1","Igf1","Csf1","Plau","Olr1","Lgals3","Cd63","Cd9","Nes","Vat1","Ccl3","Fosl1","Bach1","Ccl4","Plin2","Mrc1","Srxn1","Clic4","Rai14","Lpl","Top2a"),
                    gene.order = c("Csf1r","P2ry12","Tmem119","Itgax","Spp1","Gpnmb","Piezo1","Igf1","Csf1","Plau","Olr1","Lgals3","Cd63","Cd9","Nes","Vat1","Ccl3","Fosl1","Bach1","Ccl4","Plin2","Mrc1","Srxn1","Clic4","Rai14","Lpl","Top2a"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot.pdf", plot = jjplot, device = 'pdf', width = 40, height = 40, units = 'cm')

#---------Removal of extraneous cells
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("0","1","2","3")]
save(AB,file = "AB.RData")
#
AB <- myfunction1(AB)
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
save(AB,file="AB2.RData")
#
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:11,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:11)
AB <- FindClusters(AB, resolution = 0.4) 
save(AB,file="AB2_new.RData")

###--
p7 <- DimPlot(AB, reduction = "umap", group.by = "group", pt.size=1,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=1, label = TRUE,repel = TRUE,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2_Cluster3.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

jjplot <- jjDotPlot(object = AB,
                    gene = c("Csf1r","P2ry12","Tmem119","Itgax","Spp1","Gpnmb","Piezo1","Igf1","Csf1","Plau","Olr1","Lgals3","Cd63","Cd9","Nes","Vat1","Ccl3","Fosl1","Bach1","Ccl4","Plin2","Mrc1","Srxn1","Clic4","Rai14","Lpl","Top2a"),
                    gene.order = c("Csf1r","P2ry12","Tmem119","Itgax","Spp1","Gpnmb","Piezo1","Igf1","Csf1","Plau","Olr1","Lgals3","Cd63","Cd9","Nes","Vat1","Ccl3","Fosl1","Bach1","Ccl4","Plin2","Mrc1","Srxn1","Clic4","Rai14","Lpl","Top2a"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot3.pdf", plot = jjplot, device = 'pdf', width = 25, height = 25, units = 'cm')

table(AB@meta.data[["seurat_clusters"]], AB@meta.data[["orig.ident"]])

###--rename cluster
current.cluster.ids <- c("0","1","2","3")
new.cluster.ids <- c("SPP1L","SPP1L","SPP1L","SPP1H")
AB@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(AB@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(AB@meta.data$celltype)
table(AB@meta.data[["celltype"]], AB@meta.data[["age"]])
Idents(AB) <- factor(AB$celltype, levels = c("SPP1H","SPP1L"))
DimPlot(AB, reduction = "umap", label = TRUE)
save(AB,file="MG(SPP1L-SPP1H).RData")

#-----------MAST+vol
library(Seurat)
library(ggVolcano)
###
##MAST
library(MAST);library(Seurat);library(dplyr)
#
Idents(AB)="celltype"
table(Idents(AB)) 
#AB <- AB[,AB@meta.data[["celltype"]] %in% c("SPP1H","SPP1L")]
#AB = subset(AB,downsample=1000)
fData = data.frame(symbolid=rownames(AB),primerid=rownames(AB))
rownames(fData)=fData$symbolid
cData = AB@meta.data
cData$wellKey <- rownames(cData)
sca = FromMatrix(as.matrix(AB@assays$RNA@data), cData = cData,fData = fData)
gc()
dim(sca)
table(colData(sca)$celltype)
cond<-factor(colData(sca)$celltype)
cond<-relevel(cond,"SPP1L")
colData(sca)$condition<-cond
## (1) Calibration of cngeneson covariates: default parameters
zlmCond <- zlm(~condition + nCount_RNA + percent.mt + group , sca, method="bayesglm", ebayes=TRUE)
summaryCond <- summary(zlmCond,doLRT='conditionSPP1H')
summaryDt <- summaryCond$datatable
levels(summaryDt$contrast)

df_pval = summaryDt %>% 
  dplyr::filter(contrast=='conditionSPP1H') %>% 
  dplyr::filter(component=='H') %>% 
  dplyr::select(primerid, `Pr(>Chisq)`)

df_logfc = summaryDt %>% 
  dplyr::filter(contrast=='conditionSPP1H') %>% 
  dplyr::filter(component=='logFC') %>% 
  dplyr::select(primerid, coef, ci.hi, ci.lo)

df_stat = dplyr::inner_join(df_logfc, df_pval) %>% 
  dplyr::rename("symbol"="primerid") %>% 
  dplyr::rename("pval"="Pr(>Chisq)","logFC"="coef") %>% 
  dplyr::mutate("fdr" = p.adjust(pval)) %>% 
  dplyr::arrange(fdr)
head(df_stat)

df_stat$FC<-10^(abs(df_stat$logFC))
df_stat$FC<-ifelse(df_stat$logFC>0,df_stat$FC*(1),df_stat$FC*-1) 
df_stat$log2FC <- log2(abs(df_stat$FC))
df_stat$log2FC <- ifelse(df_stat$FC>0,df_stat$log2FC*(1),df_stat$log2FC*-1)
write.csv(df_stat,"MAST_DEGs_SPP1H_SPP1L.csv")

markers <- read.csv2("MAST_DEGs_SPP1H_SPP1L.csv",sep=",",row.names=1)
row.names(markers) <- markers[,1]
colnames(markers)[1] <-"gene"
colnames(markers)[8] <-"log2FoldChange"
colnames(markers)[6] <-"padj"
markers$log2FoldChange <- as.numeric(markers$log2FoldChange)
markers$padj <- as.numeric(markers$padj)
#
markers <- na.omit(markers)
data <- add_regulate(markers, log2FC_name = "log2FoldChange",
                     fdr_name = "padj",log2FC = 0.25, fdr = 0.05)
# plot
ggvolcano(data, x = "log2FoldChange", y = "padj", 
          label = "gene", label_number = 20, output = FALSE)
library(RColorBrewer)
p1 <- ggvolcano(data, x = "log2FoldChange", y = "padj",
                fills = c("#607999","#E6E6DC","#B0765E"),
                colors = c("#607999","#E6E6DC","#B0765E"),
                log2FC_cut = 0.25,
                FDR_cut = 0.05,
                label = "gene", label_number = 20, custom_label = c("Myo1e","Spp1",	"Rbpj",	"Abr",	"Igf1",	"Qk", "Slc9a9", "Inpp5d"), pointSize = 3, pointShape = 17, output = FALSE)
p1
ggsave(filename = "Vol-MAST(SPP1H_SPP1L).pdf", plot = p1, device = 'pdf', width = 15, height = 15, units = 'cm')

#---------Convert to h5ad
library(Seurat)
#remotes::install_github("satijalab/seurat-data")
library(SeuratData)
library(SeuratDisk)
#InstallData("pbmc3k")
# An object of class Seurat 
# 13714 features across 2638 samples within 1 assay 
# Active assay: RNA (13714 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap
SaveH5Seurat(AB, filename = "AB.h5ad.h5seurat")
Convert("AB.h5ad.h5seurat", dest = "h5ad")

AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("3")]
save(AB,file="AB(SPP1).RData")

#-------------GSE142445-------------#
#
Mm_Sham1 <- Read10X(data.dir = "Data/Sham1/")
Mm_Sham2 <- Read10X(data.dir = "Data/Sham2/")
Mm_Sham3 <- Read10X(data.dir = "Data/Sham3/")
Mm_S4h <- Read10X(data.dir = "Data/D4h/")
Mm_S1d <- Read10X(data.dir = "Data/D1d/")
Mm_S3d <- Read10X(data.dir = "Data/D3d/")
Mm_S7d <- Read10X(data.dir = "Data/D7d/")
#
Mm_Sham1 <- CreateSeuratObject(counts = Mm_Sham1, project = "Mm_Sham1", min.cells = 3, min.features = 200)
Mm_Sham2 <- CreateSeuratObject(counts = Mm_Sham2, project = "Mm_Sham2", min.cells = 3, min.features = 200)
Mm_Sham3 <- CreateSeuratObject(counts = Mm_Sham3, project = "Mm_Sham3", min.cells = 3, min.features = 200)
Mm_S4h <- CreateSeuratObject(counts = Mm_S4h, project = "Mm_S4h", min.cells = 3, min.features = 200)
Mm_S1d <- CreateSeuratObject(counts = Mm_S1d, project = "Mm_S1d", min.cells = 3, min.features = 200)
Mm_S3d <- CreateSeuratObject(counts = Mm_S3d, project = "Mm_S3d", min.cells = 3, min.features = 200)
Mm_S7d <- CreateSeuratObject(counts = Mm_S7d, project = "Mm_S7d", min.cells = 3, min.features = 200)
#
Mm_Sham1[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham1, pattern = "^Mt-")
Mm_Sham2[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham2, pattern = "^Mt-")
Mm_Sham3[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham3, pattern = "^Mt-")
Mm_S4h[["percent.mt"]] <- PercentageFeatureSet(Mm_S4h, pattern = "^Mt-")
Mm_S1d[["percent.mt"]] <- PercentageFeatureSet(Mm_S1d, pattern = "^Mt-")
Mm_S3d[["percent.mt"]] <- PercentageFeatureSet(Mm_S3d, pattern = "^Mt-")
Mm_S7d[["percent.mt"]] <- PercentageFeatureSet(Mm_S7d, pattern = "^Mt-")
#
Mm_Sham1 <- subset(Mm_Sham1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
Mm_Sham2 <- subset(Mm_Sham2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
Mm_Sham3 <- subset(Mm_Sham3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
Mm_S4h <- subset(Mm_S4h, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
Mm_S1d <- subset(Mm_S1d, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
Mm_S3d <- subset(Mm_S3d, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
Mm_S7d <- subset(Mm_S7d, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
#
myfunction1 <- function(xxxx){
  xxxx <- NormalizeData(xxxx, normalization.method = "LogNormalize", scale.factor = 10000)
  xxxx <- FindVariableFeatures(xxxx, selection.method = "vst", nfeatures = 2000)
  return(xxxx)   }
AB <- merge(Mm_Sham1,y=c(Mm_Sham2,Mm_Sham3,Mm_S4h,Mm_S1d,Mm_S3d,Mm_S7d))
AB <- myfunction1(AB)
#
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
ggsave(filename = "ElbowPlot.pdf", plot = plot2, device = 'pdf', width = 14, height = 12, units = 'cm')
rm('plot1','plot2')
save(AB,file="AB.RData")
#
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:15,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:15)
AB <- FindClusters(AB, resolution = 0.4) 
AB$age=str_replace(AB$orig.ident,"Mm_","")
#
AB@meta.data[["age"]]<-factor(AB@meta.data[["age"]], levels=c("Sham","S1d","S3d","S7d"))
#
DefaultAssay(AB) <- "RNA"
save(AB,file="AB_new.RData")
#--
p7 <- DimPlot(AB, reduction = "umap", group.by = "orig.ident", pt.size=1,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=1, label = TRUE,repel = TRUE,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2_Cluster2.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

AB100=subset(AB,downsample=100)
markers100 <- FindAllMarkers(AB100, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox")
write.table(markers100,file="markers100.txt",quote=F,sep="\t",row.names=F,col.names=T)

#Extraction of microglia
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("5",	"11","14")]
save(AB,file="AB(MG).RData")
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
ggsave(filename = "ElbowPlot2.pdf", plot = plot2, device = 'pdf', width = 14, height = 12, units = 'cm')
rm('plot1','plot2')
save(AB,file="MG.RData")
#
#AB <- FindNeighbors(AB, dims = 1:30)
#AB <- FindClusters(AB, resolution = 0.8)
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:8,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:8)
AB <- FindClusters(AB, resolution = 0.3) 

DimPlot(AB, reduction = "umap")
save(AB,file="MG_new.RData")
#
AB$group=str_replace(AB$orig.ident,"Mm_","")
AB$group=str_replace(AB$group,"Sham.*","Sham")
AB@meta.data[["group"]]<-factor(AB@meta.data[["group"]], levels=c("Sham","S4h","S1d","S3d","S7d"))
#--
p7 <- DimPlot(AB, cols = c('#f7fc7d','#6ac626','#28bd6a','#006eb8','#090ba2'), reduction = "umap", group.by = "group", pt.size=1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, cols = c('#F0B371',"#ee695b",'#C44270','#e4c2ce','#8c6fa3','#b96f94','#C07AD2'), reduction = "umap", group.by = "ident", pt.size=1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2(MG).pdf", plot = umap2, device = 'pdf', width = 33, height = 13, units = 'cm')
#---
jjplot <- jjDotPlot(object = AB,
                    gene = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1",'Tmem119',"Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Igf1","Spp1","Ccl4","Itgb3","Csf1","Gpnmb","Plin2","Mrc1","C1qa","C1qb","C1qc","Srxn1"),
                    gene.order = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1",'Tmem119',"Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Igf1","Spp1","Ccl4","Itgb3","Csf1","Gpnmb","Plin2","Mrc1","C1qa","C1qb","C1qc","Srxn1"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
jjplot <- jjDotPlot(object = AB,
                    gene = c("C1qa","C1qb","C1qc","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1",'Tmem119',"Igf1","Spp1","Ccl4","Itgax","Csf1","Gpnmb","Plin2","Mrc1","Srxn1","Olr1","Clic4","Rai14","Ccl3","Lgals3","Lpl","Cd63"),
                    gene.order = c("C1qa","C1qb","C1qc","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1",'Tmem119',"Igf1","Spp1","Ccl4","Itgax","Csf1","Gpnmb","Plin2","Mrc1","Srxn1","Olr1","Clic4","Rai14","Ccl3","Lgals3","Lpl","Cd63"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot.pdf", plot = jjplot, device = 'pdf', width = 40, height = 40, units = 'cm')
table(AB@meta.data[["seurat_clusters"]], AB@meta.data[["orig.ident"]])
#
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("3")]
save(AB,file="AB(SPP1).RData")

library(Seurat)
#remotes::install_github("satijalab/seurat-data")
library(SeuratData)
library(SeuratDisk)
#InstallData("pbmc3k")
# An object of class Seurat 
# 13714 features across 2638 samples within 1 assay 
# Active assay: RNA (13714 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap
SaveH5Seurat(AB, filename = "AB.h5ad.h5seurat")
Convert("AB.h5ad.h5seurat", dest = "h5ad")

#-----------MAST+vol
library(Seurat)
library(ggVolcano)
###
##MAST
library(MAST);library(Seurat);library(dplyr)
#
Idents(AB)="celltype"
table(Idents(AB)) 
#AB <- AB[,AB@meta.data[["celltype"]] %in% c("SPP1H","SPP1L")]
#AB = subset(AB,downsample=1000)
fData = data.frame(symbolid=rownames(AB),primerid=rownames(AB))
rownames(fData)=fData$symbolid
cData = AB@meta.data
cData$wellKey <- rownames(cData)
sca = FromMatrix(as.matrix(AB@assays$RNA@data), cData = cData,fData = fData)
gc()
dim(sca)
table(colData(sca)$celltype)
cond<-factor(colData(sca)$celltype)
cond<-relevel(cond,"SPP1L")
colData(sca)$condition<-cond
## (1) Calibration of cngeneson covariates: default parameters
zlmCond <- zlm(~condition + nCount_RNA + percent.mt + group , sca, method="bayesglm", ebayes=TRUE)
summaryCond <- summary(zlmCond,doLRT='conditionSPP1H')
summaryDt <- summaryCond$datatable
levels(summaryDt$contrast)

df_pval = summaryDt %>% 
  dplyr::filter(contrast=='conditionSPP1H') %>% 
  dplyr::filter(component=='H') %>% 
  dplyr::select(primerid, `Pr(>Chisq)`)

df_logfc = summaryDt %>% 
  dplyr::filter(contrast=='conditionSPP1H') %>% 
  dplyr::filter(component=='logFC') %>% 
  dplyr::select(primerid, coef, ci.hi, ci.lo)

df_stat = dplyr::inner_join(df_logfc, df_pval) %>% 
  dplyr::rename("symbol"="primerid") %>% 
  dplyr::rename("pval"="Pr(>Chisq)","logFC"="coef") %>% 
  dplyr::mutate("fdr" = p.adjust(pval)) %>% 
  dplyr::arrange(fdr)
head(df_stat)

df_stat$FC<-10^(abs(df_stat$logFC))
df_stat$FC<-ifelse(df_stat$logFC>0,df_stat$FC*(1),df_stat$FC*-1) 
df_stat$log2FC <- log2(abs(df_stat$FC))
df_stat$log2FC <- ifelse(df_stat$FC>0,df_stat$log2FC*(1),df_stat$log2FC*-1)
write.csv(df_stat,"MAST_DEGs_SPP1H_SPP1L.csv")

markers <- read.csv2("MAST_DEGs_SPP1H_SPP1L.csv",sep=",",row.names=1)
row.names(markers) <- markers[,1]
colnames(markers)[1] <-"gene"
colnames(markers)[8] <-"log2FoldChange"
colnames(markers)[6] <-"padj"
markers$log2FoldChange <- as.numeric(markers$log2FoldChange)
markers$padj <- as.numeric(markers$padj)
#
markers <- na.omit(markers)
data <- add_regulate(markers, log2FC_name = "log2FoldChange",
                     fdr_name = "padj",log2FC = 0.25, fdr = 0.05)
# plot
ggvolcano(data, x = "log2FoldChange", y = "padj", 
          label = "gene", label_number = 20, output = FALSE)
library(RColorBrewer)
p1 <- ggvolcano(data, x = "log2FoldChange", y = "padj",
                fills = c("#607999","#E6E6DC","#B0765E"),
                colors = c("#607999","#E6E6DC","#B0765E"),
                log2FC_cut = 0.25,
                FDR_cut = 0.05,
                label = "gene", label_number = 20, custom_label = c("Ccl4",	"Ctsb",	"Lpl",	"Ftl1",	"Spp1",	"Ccl9",	"Ccl3",	"Cd14",	"Lgals3",	"Dab2",	"Cstb",	"Il1rn",	"Cd72",	"Pfn1",	"Lilrb4a",	"Cd81",	"Hsp90ab1",	"Tmem119",	"Bhlhe41",	"Ccl12"), pointSize = 3, pointShape = 17, output = FALSE)
p1
ggsave(filename = "Vol-MAST(SPP1H_SPP1L).pdf", plot = p1, device = 'pdf', width = 15, height = 15, units = 'cm')

jjplot <- jjDotPlot(object = AB,
                    gene = c("Csf1r","P2ry12","Tmem119","Spp1","Igf1","Csf1","Lgals3","Cd14",	"Adam8","Nes","Vat1","Ccl4",	"Lpl",	"Ccl3",	"Ccl9",	"Lilrb4a","Lilr4b",	"Plin2","Ell2"),
                    gene.order = c("Csf1r","P2ry12","Tmem119","Spp1","Igf1","Csf1","Lgals3","Cd14",	"Adam8","Nes","Vat1","Ccl4",	"Lpl",	"Ccl3",	"Ccl9",	"Lilrb4a","Lilr4b",	"Plin2","Ell2"),
                    id = 'celltype',
                    cluster.order = c("0",	"1",	"2",	"3",	"4",	"5"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot-cluster.pdf", plot = jjplot, device = 'pdf', width = 22, height = 11, units = 'cm')

#-------------GSE174574-------------#
#
Mm_Sham1 <- Read10X(data.dir = "Data/Sham1/")
Mm_Sham2 <- Read10X(data.dir = "Data/Sham2/")
Mm_Sham3 <- Read10X(data.dir = "Data/Sham3/")
Mm_S1a <- Read10X(data.dir = "Data/D1a/")
Mm_S1b <- Read10X(data.dir = "Data/D1b/")
Mm_S1c <- Read10X(data.dir = "Data/D1c/")
#
Mm_Sham1 <- CreateSeuratObject(counts = Mm_Sham1, project = "Mm_Sham1", min.cells = 3, min.features = 200)
Mm_Sham2 <- CreateSeuratObject(counts = Mm_Sham2, project = "Mm_Sham2", min.cells = 3, min.features = 200)
Mm_Sham3 <- CreateSeuratObject(counts = Mm_Sham3, project = "Mm_Sham3", min.cells = 3, min.features = 200)
Mm_S1a <- CreateSeuratObject(counts = Mm_S1a, project = "Mm_S1a", min.cells = 3, min.features = 200)
Mm_S1b <- CreateSeuratObject(counts = Mm_S1b, project = "Mm_S1b", min.cells = 3, min.features = 200)
Mm_S1c <- CreateSeuratObject(counts = Mm_S1c, project = "Mm_S1c", min.cells = 3, min.features = 200)
#
Mm_Sham1[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham1, pattern = "^Mt-")
Mm_Sham2[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham2, pattern = "^Mt-")
Mm_Sham3[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham3, pattern = "^Mt-")
Mm_S1a[["percent.mt"]] <- PercentageFeatureSet(Mm_S1a, pattern = "^Mt-")
Mm_S1b[["percent.mt"]] <- PercentageFeatureSet(Mm_S1b, pattern = "^Mt-")
Mm_S1c[["percent.mt"]] <- PercentageFeatureSet(Mm_S1c, pattern = "^Mt-")
#
Mm_Sham1 <- subset(Mm_Sham1, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
Mm_Sham2 <- subset(Mm_Sham2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
Mm_Sham3 <- subset(Mm_Sham3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
Mm_S1a <- subset(Mm_S1a, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
Mm_S1b <- subset(Mm_S1b, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)
Mm_S1c <- subset(Mm_S1c, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 25)

myfunction1 <- function(xxxx){
  xxxx <- NormalizeData(xxxx, normalization.method = "LogNormalize", scale.factor = 10000)
  xxxx <- FindVariableFeatures(xxxx, selection.method = "vst", nfeatures = 2000)
  return(xxxx)   }
AB <- merge(Mm_Sham1,y=c(Mm_Sham2,Mm_Sham3,Mm_S1a,Mm_S1b,Mm_S1c))
AB <- myfunction1(AB)
#
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
ggsave(filename = "ElbowPlot.pdf", plot = plot2, device = 'pdf', width = 14, height = 12, units = 'cm')
rm('plot1','plot2')
save(AB,file="AB.RData")
#
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:30,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:30)
AB <- FindClusters(AB, resolution = 0.4) 
AB$age=str_replace(AB$orig.ident,"Mm_","")
save(AB,file="AB_new.RData")
#--
p7 <- DimPlot(AB, reduction = "umap", group.by = "orig.ident", pt.size=0.01,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=0.1, label = TRUE,repel = TRUE,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2_Cluster2.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')
#
AB100=subset(AB,downsample=100)
markers100 <- FindAllMarkers(AB100, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox") 
write.table(markers100,file="markers100.txt",quote=F,sep="\t",row.names=F,col.names=T)
#
jjplot <- jjDotPlot(object = AB,
                    gene = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1","Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Igf1","Spp1","Ccl4","Itgb3","Plin2","Csf1","Gpnmb","Mrc1","C1qa","C1qb","C1qc","Itgb2","Srxn1"),
                    gene.order = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1","Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Igf1","Spp1","Ccl4","Itgb3","Plin2","Csf1","Gpnmb","Mrc1","C1qa","C1qb","C1qc","Itgb2","Srxn1"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot.pdf", plot = jjplot, device = 'pdf', width = 40, height = 40, units = 'cm')

#--------------Extraction of microglia
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("1", "3")]
save(AB,file="AB(MG).RData")
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
ggsave(filename = "ElbowPlot2.pdf", plot = plot2, device = 'pdf', width = 14, height = 12, units = 'cm')
rm('plot1','plot2')
save(AB,file="MG.RData")
#
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:11,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:11)
AB <- FindClusters(AB, resolution = 0.4) 
#AB <- RunUMAP(AB, dims = 1:30)
DimPlot(AB, reduction = "umap")
#-------Removing stray cells
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("0", "1","2","3","4","5")]
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
plot2
save(AB,file="MG2.RData")
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:5,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:5)
AB <- FindClusters(AB, resolution = 0.4) 
DimPlot(AB, reduction = "umap")
save(AB,file="MG_new.RData")
#
AB$group=str_replace(AB$orig.ident,"Mm_","")
AB$group=str_replace(AB$group,"Sham.*","Sham")
AB$group=str_replace(AB$group,"S1.*","S1d")
#
AB@meta.data[["group"]]<-factor(AB@meta.data[["group"]], levels=c("Sham","S1d"))
###--
p7 <- DimPlot(AB, cols = c('#f7fc7d','#6ac626','#28bd6a','#006eb8','#090ba2'), reduction = "umap", group.by = "group", pt.size=1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, cols = c('#F0B371',"#ee695b",'#C44270','#e4c2ce','#8c6fa3','#b96f94','#C07AD2'), reduction = "umap", group.by = "ident", pt.size=1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2(MG).pdf", plot = umap2, device = 'pdf', width = 33, height = 13, units = 'cm')
#
jjplot <- jjDotPlot(object = AB,
                    gene = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1",'Tmem119',"Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Igf1","Spp1","Ccl4","Itgb3","Csf1","Gpnmb","Plin2","Mrc1","C1qa","C1qb","C1qc","Srxn1"),
                    gene.order = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1",'Tmem119',"Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Igf1","Spp1","Ccl4","Itgb3","Csf1","Gpnmb","Plin2","Mrc1","C1qa","C1qb","C1qc","Srxn1"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
#
jjplot <- jjDotPlot(object = AB,
                    gene = c("C1qa","C1qb","C1qc","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1",'Tmem119',"Igf1","Spp1","Ccl4","Itgax","Csf1","Gpnmb","Plin2","Mrc1","Srxn1","Olr1","Clic4","Rai14","Ccl3","Lgals3","Lpl","Cd63"),
                    gene.order = c("C1qa","C1qb","C1qc","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1",'Tmem119',"Igf1","Spp1","Ccl4","Itgax","Csf1","Gpnmb","Plin2","Mrc1","Srxn1","Olr1","Clic4","Rai14","Ccl3","Lgals3","Lpl","Cd63"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot.pdf", plot = jjplot, device = 'pdf', width = 40, height = 40, units = 'cm')
#
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("2")]
save(AB,file="AB(SPP1).RData")
#
library(Seurat)
#remotes::install_github("satijalab/seurat-data")
library(SeuratData)
library(SeuratDisk)
#InstallData("pbmc3k")
# An object of class Seurat 
# 13714 features across 2638 samples within 1 assay 
# Active assay: RNA (13714 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap
SaveH5Seurat(AB, filename = "AB.h5ad.h5seurat")
Convert("AB.h5ad.h5seurat", dest = "h5ad")
#
jjplot <- jjDotPlot(object = AB,
                    gene = c("Csf1r","P2ry12","Tmem119","Spp1","Igf1","Csf1","Lgals3","Cd14",	"Mmp12","Il1rn","Cxcl2","Adam8","Vat1",	"Lpl",	"Ccl9",	"Lilrb4a","Lilr4b",	"Plin2","Ell2"),
                    gene.order = c("Csf1r","P2ry12","Tmem119","Spp1","Igf1","Csf1","Lgals3","Cd14",	"Mmp12","Il1rn","Cxcl2","Adam8","Vat1",	"Lpl",	"Ccl9",	"Lilrb4a","Lilr4b",	"Plin2","Ell2"),
                    id = 'celltype',
                    cluster.order = c("0",	"1",	"2",	"3",	"4",	"5"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)#+ggplot2:::coord_flip()
jjplot

ggsave(filename = "jjDotPlot-cluster.pdf", plot = jjplot, device = 'pdf', width = 22, height = 11, units = 'cm')

#-----------MAST+vol
library(Seurat)
library(ggVolcano)
###
##MAST
library(MAST);library(Seurat);library(dplyr)
#
Idents(AB)="celltype"
table(Idents(AB)) 
#AB <- AB[,AB@meta.data[["celltype"]] %in% c("SPP1H","SPP1L")]
#AB = subset(AB,downsample=1000)
fData = data.frame(symbolid=rownames(AB),primerid=rownames(AB))
rownames(fData)=fData$symbolid
cData = AB@meta.data
cData$wellKey <- rownames(cData)
sca = FromMatrix(as.matrix(AB@assays$RNA@data), cData = cData,fData = fData)
gc()
dim(sca)
table(colData(sca)$celltype)
cond<-factor(colData(sca)$celltype)
cond<-relevel(cond,"SPP1L")
colData(sca)$condition<-cond
## (1) Calibration of cngeneson covariates: default parameters
zlmCond <- zlm(~condition + nCount_RNA + percent.mt + group , sca, method="bayesglm", ebayes=TRUE)
summaryCond <- summary(zlmCond,doLRT='conditionSPP1H')
summaryDt <- summaryCond$datatable
levels(summaryDt$contrast)

df_pval = summaryDt %>% 
  dplyr::filter(contrast=='conditionSPP1H') %>% 
  dplyr::filter(component=='H') %>% 
  dplyr::select(primerid, `Pr(>Chisq)`)

df_logfc = summaryDt %>% 
  dplyr::filter(contrast=='conditionSPP1H') %>% 
  dplyr::filter(component=='logFC') %>% 
  dplyr::select(primerid, coef, ci.hi, ci.lo)

df_stat = dplyr::inner_join(df_logfc, df_pval) %>% 
  dplyr::rename("symbol"="primerid") %>% 
  dplyr::rename("pval"="Pr(>Chisq)","logFC"="coef") %>% 
  dplyr::mutate("fdr" = p.adjust(pval)) %>% 
  dplyr::arrange(fdr)
head(df_stat)

df_stat$FC<-10^(abs(df_stat$logFC))
df_stat$FC<-ifelse(df_stat$logFC>0,df_stat$FC*(1),df_stat$FC*-1) 
df_stat$log2FC <- log2(abs(df_stat$FC))
df_stat$log2FC <- ifelse(df_stat$FC>0,df_stat$log2FC*(1),df_stat$log2FC*-1)
write.csv(df_stat,"MAST_DEGs_SPP1H_SPP1L.csv")

markers <- read.csv2("MAST_DEGs_SPP1H_SPP1L.csv",sep=",",row.names=1)
row.names(markers) <- markers[,1]
colnames(markers)[1] <-"gene"
colnames(markers)[8] <-"log2FoldChange"
colnames(markers)[6] <-"padj"
markers$log2FoldChange <- as.numeric(markers$log2FoldChange)
markers$padj <- as.numeric(markers$padj)
#
markers <- na.omit(markers)
data <- add_regulate(markers, log2FC_name = "log2FoldChange",
                     fdr_name = "padj",log2FC = 0.25, fdr = 0.05)
# plot
ggvolcano(data, x = "log2FoldChange", y = "padj", 
          label = "gene", label_number = 20, output = FALSE)
library(RColorBrewer)
p1 <- ggvolcano(data, x = "log2FoldChange", y = "padj",
                fills = c("#607999","#E6E6DC","#B0765E"),
                colors = c("#607999","#E6E6DC","#B0765E"),
                log2FC_cut = 0.25,
                FDR_cut = 0.05,
                label = "gene", label_number = 20, custom_label = c("Lilrb4a",	"Lgals3",	"Cd14",	"Lilr4b",	"Spp1","Tmem119",	"Selplg",	"Cd81"), pointSize = 3, pointShape = 17, output = FALSE)
p1
ggsave(filename = "Vol-MAST(SPP1H_SPP1L).pdf", plot = p1, device = 'pdf', width = 15, height = 15, units = 'cm')

#-------------GSE225948-------------#
#
Mm_Sham1 <- read.csv("data/Sham1/GSM7060815_Brain_GR180716_counts.csv", sep=",", header = T,row.names = 1)
Mm_Sham2 <- read.csv("data/Sham2/GSM7060816_Brain_GR181128_counts.csv", sep=",", header = T,row.names = 1)
Mm_Sham3 <- read.csv("data/Sham3/GSM7060817_Brain_GR181212_counts.csv", sep=",", header = T,row.names = 1)
Mm_Sham4 <- read.csv("data/Sham4/GSM7060818_Brain_GR190110_counts.csv", sep=",", header = T,row.names = 1)
Mm_Mcao2d1 <- read.csv("data/Mcao2d-1/GSM7060819_Brain_GR180426_counts.csv", sep=",", header = T,row.names = 1)
Mm_Mcao2d2 <- read.csv("data/Mcao2d-2/GSM7060820_Brain_GR180614_counts.csv", sep=",", header = T,row.names = 1)
Mm_Mcao2d3 <- read.csv("data/Mcao2d-3/GSM7060821_Brain_GR180919_counts.csv", sep=",", header = T,row.names = 1)
Mm_Mcao2d4 <- read.csv("data/Mcao2d-4/GSM7060822_Brain_GR181024_counts.csv", sep=",", header = T,row.names = 1)
Mm_Mcao14d1 <- read.csv("data/Mcao14d-1/GSM7060823_Brain_GR180125_counts.csv", sep=",", header = T,row.names = 1)
Mm_Mcao14d2 <- read.csv("data/Mcao14d-2/GSM7060824_Brain_GR180613_counts.csv", sep=",", header = T,row.names = 1)
Mm_Mcao14d3 <- read.csv("data/Mcao14d-3/GSM7060825_Brain_GR180905_counts.csv", sep=",", header = T,row.names = 1)
Mm_Mcao14d4 <- read.csv("data/Mcao14d-4/GSM7060826_Brain_GR181114_counts.csv", sep=",", header = T,row.names = 1)
#
Mm_Sham1 <- CreateSeuratObject(counts = Mm_Sham1, project = "Mm_Sham1", min.cells = 3, min.features = 200)
Mm_Sham2 <- CreateSeuratObject(counts = Mm_Sham2, project = "Mm_Sham2", min.cells = 3, min.features = 200)
Mm_Sham3 <- CreateSeuratObject(counts = Mm_Sham3, project = "Mm_Sham3", min.cells = 3, min.features = 200)
Mm_Sham4 <- CreateSeuratObject(counts = Mm_Sham4, project = "Mm_Sham4", min.cells = 3, min.features = 200)
Mm_Mcao2d1 <- CreateSeuratObject(counts = Mm_Mcao2d1, project = "Mm_Mcao2d1", min.cells = 3, min.features = 200)
Mm_Mcao2d2 <- CreateSeuratObject(counts = Mm_Mcao2d2, project = "Mm_Mcao2d2", min.cells = 3, min.features = 200)
Mm_Mcao2d3 <- CreateSeuratObject(counts = Mm_Mcao2d3, project = "Mm_Mcao2d3", min.cells = 3, min.features = 200)
Mm_Mcao2d4 <- CreateSeuratObject(counts = Mm_Mcao2d4, project = "Mm_Mcao2d4", min.cells = 3, min.features = 200)
Mm_Mcao14d1 <- CreateSeuratObject(counts = Mm_Mcao14d1, project = "Mm_Mcao14d1", min.cells = 3, min.features = 200)
Mm_Mcao14d2 <- CreateSeuratObject(counts = Mm_Mcao14d2, project = "Mm_Mcao14d2", min.cells = 3, min.features = 200)
Mm_Mcao14d3 <- CreateSeuratObject(counts = Mm_Mcao14d3, project = "Mm_Mcao14d3", min.cells = 3, min.features = 200)
Mm_Mcao14d4 <- CreateSeuratObject(counts = Mm_Mcao14d4, project = "Mm_Mcao14d4", min.cells = 3, min.features = 200)
#
Mm_Sham1[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham1, pattern = "^Mt-")
Mm_Sham2[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham2, pattern = "^Mt-")
Mm_Sham3[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham3, pattern = "^Mt-")
Mm_Sham4[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham4, pattern = "^Mt-")
Mm_Mcao2d1[["percent.mt"]] <- PercentageFeatureSet(Mm_Mcao2d1, pattern = "^Mt-")
Mm_Mcao2d2[["percent.mt"]] <- PercentageFeatureSet(Mm_Mcao2d2, pattern = "^Mt-")
Mm_Mcao2d3[["percent.mt"]] <- PercentageFeatureSet(Mm_Mcao2d3, pattern = "^Mt-")
Mm_Mcao2d4[["percent.mt"]] <- PercentageFeatureSet(Mm_Mcao2d4, pattern = "^Mt-")
Mm_Mcao14d1[["percent.mt"]] <- PercentageFeatureSet(Mm_Mcao14d1, pattern = "^Mt-")
Mm_Mcao14d2[["percent.mt"]] <- PercentageFeatureSet(Mm_Mcao14d2, pattern = "^Mt-")
Mm_Mcao14d3[["percent.mt"]] <- PercentageFeatureSet(Mm_Mcao14d3, pattern = "^Mt-")
Mm_Mcao14d4[["percent.mt"]] <- PercentageFeatureSet(Mm_Mcao14d4, pattern = "^Mt-")
#
Mm_Sham1 <- subset(Mm_Sham1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Sham2 <- subset(Mm_Sham2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Sham3 <- subset(Mm_Sham3, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Sham4 <- subset(Mm_Sham4, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Mcao2d1 <- subset(Mm_Mcao2d1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Mcao2d2 <- subset(Mm_Mcao2d2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Mcao2d3 <- subset(Mm_Mcao2d3, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Mcao2d4 <- subset(Mm_Mcao2d4, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Mcao14d1 <- subset(Mm_Mcao14d1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Mcao14d2 <- subset(Mm_Mcao14d2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Mcao14d3 <- subset(Mm_Mcao14d3, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
Mm_Mcao14d4 <- subset(Mm_Mcao14d4, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
#
myfunction1 <- function(xxxx){
  xxxx <- NormalizeData(xxxx, normalization.method = "LogNormalize", scale.factor = 10000)
  xxxx <- FindVariableFeatures(xxxx, selection.method = "vst", nfeatures = 2000)
  return(xxxx)   }
AB <- merge(Mm_Sham1,y=c(Mm_Sham2,Mm_Sham3,Mm_Sham4,Mm_Mcao2d1,Mm_Mcao2d2,Mm_Mcao2d3,Mm_Mcao2d4,Mm_Mcao14d1,Mm_Mcao14d2,Mm_Mcao14d3,Mm_Mcao14d4))
AB <- myfunction1(AB)
#
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
ggsave(filename = "ElbowPlot.pdf", plot = plot2, device = 'pdf', width = 14, height = 12, units = 'cm')
rm('plot1','plot2')
save(AB,file="AB.RData")
#
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:30,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:30)
AB <- FindClusters(AB, resolution = 0.8) 
AB$age=str_replace(AB$orig.ident,"Mm_","")
#
save(AB,file="AB_new.RData")
#
p7 <- DimPlot(AB, reduction = "umap", group.by = "orig.ident", pt.size=0.01,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2_Cluster2.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

AB100=subset(AB,downsample=100)
markers100 <- FindAllMarkers(AB100, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox")  ##耗时久
write.table(markers100,file="markers100.txt",quote=F,sep="\t",row.names=F,col.names=T)
#Extraction of microglia
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("0", "1", "5", "7", "15", "18")]
save(AB,file="AB(MG).RData")
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
plot2
save(AB,file="MG.RData")
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:10,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:10)
AB <- FindClusters(AB, resolution = 0.3) 
#
AB$group=str_replace(AB$orig.ident,"Mm_","")
AB$group=str_replace(AB$group,"Sham.*","Sham")
AB$group=str_replace(AB$group,"Mcao","S")
AB$group=str_replace(AB$group,"S14d.*","S14d")
AB$group=str_replace(AB$group,"S2d.*","S2d")
save(AB,file="MG_new.RData")
###--
p7 <- DimPlot(AB, cols = c('#f7fc7d','#28bd6a','#090ba2'), reduction = "umap", group.by = "group", pt.size=1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, cols = c('#F0B371',"#ee695b",'#C44270','#e4c2ce','#8c6fa3','#b96f94','#C07AD2'), reduction = "umap", group.by = "ident", pt.size=1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2(MG).pdf", plot = umap2, device = 'pdf', width = 33, height = 13, units = 'cm')
#
library(Seurat)
#remotes::install_github("satijalab/seurat-data")
library(SeuratData)
library(SeuratDisk)
#InstallData("pbmc3k")
# An object of class Seurat 
# 13714 features across 2638 samples within 1 assay 
# Active assay: RNA (13714 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap
SaveH5Seurat(AB, filename = "AB.h5ad.h5seurat")
Convert("AB.h5ad.h5seurat", dest = "h5ad")
#choice SPP1
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("2")]
save(AB,file="AB(SPP1).RData")

#-----------MAST+vol
library(Seurat)
library(ggVolcano)
###
##MAST
library(MAST);library(Seurat);library(dplyr)
#
Idents(AB)="celltype"
table(Idents(AB)) 
#AB <- AB[,AB@meta.data[["celltype"]] %in% c("SPP1H","SPP1L")]
#AB = subset(AB,downsample=1000)
fData = data.frame(symbolid=rownames(AB),primerid=rownames(AB))
rownames(fData)=fData$symbolid
cData = AB@meta.data
cData$wellKey <- rownames(cData)
sca = FromMatrix(as.matrix(AB@assays$RNA@data), cData = cData,fData = fData)
gc()
dim(sca)
table(colData(sca)$celltype)
cond<-factor(colData(sca)$celltype)
cond<-relevel(cond,"SPP1L")
colData(sca)$condition<-cond
## (1) Calibration of cngeneson covariates: default parameters
zlmCond <- zlm(~condition + nCount_RNA + percent.mt + group , sca, method="bayesglm", ebayes=TRUE)
summaryCond <- summary(zlmCond,doLRT='conditionSPP1H')
summaryDt <- summaryCond$datatable
levels(summaryDt$contrast)

df_pval = summaryDt %>% 
  dplyr::filter(contrast=='conditionSPP1H') %>% 
  dplyr::filter(component=='H') %>% 
  dplyr::select(primerid, `Pr(>Chisq)`)

df_logfc = summaryDt %>% 
  dplyr::filter(contrast=='conditionSPP1H') %>% 
  dplyr::filter(component=='logFC') %>% 
  dplyr::select(primerid, coef, ci.hi, ci.lo)

df_stat = dplyr::inner_join(df_logfc, df_pval) %>% 
  dplyr::rename("symbol"="primerid") %>% 
  dplyr::rename("pval"="Pr(>Chisq)","logFC"="coef") %>% 
  dplyr::mutate("fdr" = p.adjust(pval)) %>% 
  dplyr::arrange(fdr)
head(df_stat)

df_stat$FC<-10^(abs(df_stat$logFC))
df_stat$FC<-ifelse(df_stat$logFC>0,df_stat$FC*(1),df_stat$FC*-1) 
df_stat$log2FC <- log2(abs(df_stat$FC))
df_stat$log2FC <- ifelse(df_stat$FC>0,df_stat$log2FC*(1),df_stat$log2FC*-1)
write.csv(df_stat,"MAST_DEGs_SPP1H_SPP1L.csv")

markers <- read.csv2("MAST_DEGs_SPP1H_SPP1L.csv",sep=",",row.names=1)
row.names(markers) <- markers[,1]
colnames(markers)[1] <-"gene"
colnames(markers)[8] <-"log2FoldChange"
colnames(markers)[6] <-"padj"
markers$log2FoldChange <- as.numeric(markers$log2FoldChange)
markers$padj <- as.numeric(markers$padj)
#
markers <- na.omit(markers)
data <- add_regulate(markers, log2FC_name = "log2FoldChange",
                     fdr_name = "padj",log2FC = 0.25, fdr = 0.05)
# plot
ggvolcano(data, x = "log2FoldChange", y = "padj", 
          label = "gene", label_number = 20, output = FALSE)
library(RColorBrewer)
p1 <- ggvolcano(data, x = "log2FoldChange", y = "padj",
                fills = c("#607999","#E6E6DC","#B0765E"),
                colors = c("#607999","#E6E6DC","#B0765E"),
                log2FC_cut = 0.25,
                FDR_cut = 0.05,
                label = "gene", label_number = 20, custom_label = c("Clec7a",	"Apoe",	"Fth1",	"Ctsb",	"Lyz2",	"Siglech",	"Cx3cr1",	"P2ry12"), pointSize = 3, pointShape = 17, output = FALSE)
p1
ggsave(filename = "Vol-MAST(SPP1H_SPP1L).pdf", plot = p1, device = 'pdf', width = 15, height = 15, units = 'cm')

#-------------SRP465344-------------#
#
Mm_Sham <- Read10X(data.dir = "Data/Sham/filtered_feature_bc_matrix/")
Mm_S3h <- Read10X(data.dir = "Data/S3h/filtered_feature_bc_matrix/")
Mm_S12h <- Read10X(data.dir = "Data/S12h/filtered_feature_bc_matrix/")
Mm_S3d <- Read10X(data.dir = "Data/S3d/filtered_feature_bc_matrix/")
#
Mm_Sham <- CreateSeuratObject(counts = Mm_Sham, project = "Mm_Sham")
Mm_S3h <- CreateSeuratObject(counts = Mm_S3h, project = "Mm_S3h")
Mm_S12h <- CreateSeuratObject(counts = Mm_S12h, project = "Mm_S12h")
Mm_S3d <- CreateSeuratObject(counts = Mm_S3d, project = "Mm_S3d")
#
Mm_Sham <- CreateSeuratObject(counts = Mm_Sham@assays$RNA@counts, project = "Mm_Sham", min.cells = 3, min.features = 200)
Mm_S3h <- CreateSeuratObject(counts = Mm_S3h@assays$RNA@counts, project = "Mm_S3h", min.cells = 3, min.features = 200)
Mm_S12h <- CreateSeuratObject(counts = Mm_S12h@assays$RNA@counts, project = "Mm_S12h", min.cells = 3, min.features = 200)
Mm_S3d <- CreateSeuratObject(counts = Mm_S3d@assays$RNA@counts, project = "Mm_S3d", min.cells = 3, min.features = 200)
#
Mm_Sham[["percent.mt"]] <- PercentageFeatureSet(Mm_Sham, pattern = "^Mt-")
Mm_S3h[["percent.mt"]] <- PercentageFeatureSet(Mm_S3h, pattern = "^Mt-")
Mm_S12h[["percent.mt"]] <- PercentageFeatureSet(Mm_S12h, pattern = "^Mt-")
Mm_S3d[["percent.mt"]] <- PercentageFeatureSet(Mm_S3d, pattern = "^Mt-")
#
Mm_Sham <- subset(Mm_Sham, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
Mm_S3h <- subset(Mm_S3h, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
Mm_S12h <- subset(Mm_S12h, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
Mm_S3d <- subset(Mm_S3d, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
#
myfunction1 <- function(xxxx){
  xxxx <- NormalizeData(xxxx, normalization.method = "LogNormalize", scale.factor = 10000)
  xxxx <- FindVariableFeatures(xxxx, selection.method = "vst", nfeatures = 2000)
  return(xxxx)   }
AB <- merge(Mm_Sham,y=c(Mm_S3h,Mm_S12h,Mm_S3d))
AB <- myfunction1(AB)
#
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
ggsave(filename = "ElbowPlot.pdf", plot = plot2, device = 'pdf', width = 14, height = 12, units = 'cm')
rm('plot1','plot2')
save(AB,file="AB.RData")
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:20,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:20)
AB <- FindClusters(AB, resolution = 0.4) 
AB@meta.data[["age"]]<-factor(AB@meta.data[["age"]], levels=c("Sham","S1d","S3d","S7d"))
save(AB,file="AB_new.RData")
###--
p7 <- DimPlot(AB, reduction = "umap", group.by = "orig.ident", pt.size=0.01,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=1, label = TRUE,repel = TRUE,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2_Cluster2.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')
#
#Extraction of microglia
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("0",	"8")]
save(AB,file="AB(MG).RData")
#
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
ggsave(filename = "ElbowPlot2.pdf", plot = plot2, device = 'pdf', width = 14, height = 12, units = 'cm')
rm('plot1','plot2')
save(AB,file="MG.RData")
#
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:20,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:20)
AB <- FindClusters(AB, resolution = 0.4) 
#
DimPlot(AB, reduction = "umap")
#
save(AB,file="MG_new.RData")
#Extraction of microglia
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("0","1","2","3","4","6","7")]
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
plot2
save(AB,file="MG2.RData")
#
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:14,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:14)
AB <- FindClusters(AB, resolution = 0.4) 
DimPlot(AB, reduction = "umap")
save(AB,file="MG_new2.RData")
#Extraction of microglia
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("0","1","2","3")]
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
plot2 <- ElbowPlot(AB)
plot2
save(AB,file="MG3.RData")
library(harmony)
AB <- RunHarmony(AB,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")
AB <- RunUMAP(AB, reduction = "harmony", dims = 1:10,reduction.name = "umap")
AB <- FindNeighbors(AB,reduction = "harmony",dims = 1:10)
AB <- FindClusters(AB, resolution = 0.4) 
DimPlot(AB, reduction = "umap")
save(AB,file="MG_new3.RData")
#
AB$group=str_replace(AB$orig.ident,"Mm_","")
###--
p7 <- DimPlot(AB, cols = c('#f7fc7d','#6ac626','#28bd6a','#006eb8','#090ba2'), reduction = "umap", group.by = "group", pt.size=1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, cols = c('#F0B371',"#ee695b",'#C44270','#e4c2ce','#8c6fa3','#b96f94','#C07AD2'), reduction = "umap", group.by = "ident", pt.size=1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2(MG).pdf", plot = umap2, device = 'pdf', width = 33, height = 13, units = 'cm')
#
library(Seurat)
#remotes::install_github("satijalab/seurat-data")
library(SeuratData)
library(SeuratDisk)
#InstallData("pbmc3k")
# An object of class Seurat 
# 13714 features across 2638 samples within 1 assay 
# Active assay: RNA (13714 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap
SaveH5Seurat(AB, filename = "AB.h5ad.h5seurat")
Convert("AB.h5ad.h5seurat", dest = "h5ad")
#
AB <- AB[,AB@meta.data[["seurat_clusters"]] %in% c("4")]
save(AB,file="AB(SPP1).RData")
#-----------MAST+vol
library(Seurat)
library(ggVolcano)
###
##MAST
library(MAST);library(Seurat);library(dplyr)
#
Idents(AB)="celltype"
table(Idents(AB)) 
#AB <- AB[,AB@meta.data[["celltype"]] %in% c("SPP1H","SPP1L")]
#AB = subset(AB,downsample=1000)
fData = data.frame(symbolid=rownames(AB),primerid=rownames(AB))
rownames(fData)=fData$symbolid
cData = AB@meta.data
cData$wellKey <- rownames(cData)
sca = FromMatrix(as.matrix(AB@assays$RNA@data), cData = cData,fData = fData)
gc()
dim(sca)
table(colData(sca)$celltype)
cond<-factor(colData(sca)$celltype)
cond<-relevel(cond,"SPP1L")
colData(sca)$condition<-cond
## (1) Calibration of cngeneson covariates: default parameters
zlmCond <- zlm(~condition + nCount_RNA + percent.mt + group , sca, method="bayesglm", ebayes=TRUE)
summaryCond <- summary(zlmCond,doLRT='conditionSPP1H')
summaryDt <- summaryCond$datatable
levels(summaryDt$contrast)

df_pval = summaryDt %>% 
  dplyr::filter(contrast=='conditionSPP1H') %>% 
  dplyr::filter(component=='H') %>% 
  dplyr::select(primerid, `Pr(>Chisq)`)

df_logfc = summaryDt %>% 
  dplyr::filter(contrast=='conditionSPP1H') %>% 
  dplyr::filter(component=='logFC') %>% 
  dplyr::select(primerid, coef, ci.hi, ci.lo)

df_stat = dplyr::inner_join(df_logfc, df_pval) %>% 
  dplyr::rename("symbol"="primerid") %>% 
  dplyr::rename("pval"="Pr(>Chisq)","logFC"="coef") %>% 
  dplyr::mutate("fdr" = p.adjust(pval)) %>% 
  dplyr::arrange(fdr)
head(df_stat)

df_stat$FC<-10^(abs(df_stat$logFC))
df_stat$FC<-ifelse(df_stat$logFC>0,df_stat$FC*(1),df_stat$FC*-1) 
df_stat$log2FC <- log2(abs(df_stat$FC))
df_stat$log2FC <- ifelse(df_stat$FC>0,df_stat$log2FC*(1),df_stat$log2FC*-1)
write.csv(df_stat,"MAST_DEGs_SPP1H_SPP1L.csv")

markers <- read.csv2("MAST_DEGs_SPP1H_SPP1L.csv",sep=",",row.names=1)
row.names(markers) <- markers[,1]
colnames(markers)[1] <-"gene"
colnames(markers)[8] <-"log2FoldChange"
colnames(markers)[6] <-"padj"
markers$log2FoldChange <- as.numeric(markers$log2FoldChange)
markers$padj <- as.numeric(markers$padj)
#
markers <- na.omit(markers)
data <- add_regulate(markers, log2FC_name = "log2FoldChange",
                     fdr_name = "padj",log2FC = 0.25, fdr = 0.05)
# plot
ggvolcano(data, x = "log2FoldChange", y = "padj", 
          label = "gene", label_number = 20, output = FALSE)
library(RColorBrewer)
p1 <- ggvolcano(data, x = "log2FoldChange", y = "padj",
                fills = c("#607999","#E6E6DC","#B0765E"),
                colors = c("#607999","#E6E6DC","#B0765E"),
                log2FC_cut = 0.25,
                FDR_cut = 0.05,
                label = "gene", label_number = 20, custom_label = c("Cd14",	"Fmnl2",	"Ccl4",	"Id2",	"Csf1",	"Selplg",	"Lgmn",	"Ptpro"), pointSize = 3, pointShape = 17, output = FALSE)
p1
ggsave(filename = "Vol-MAST(SPP1H_SPP1L).pdf", plot = p1, device = 'pdf', width = 15, height = 15, units = 'cm')

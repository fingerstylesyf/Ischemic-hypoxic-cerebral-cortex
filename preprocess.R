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

#——————————Data Integration——————————#
###
SD_Sham <- Read10X(data.dir = "Data/SD-Sham_count/")
SD_Sham3d <- Read10X(data.dir = "Data/SD-Sham3d_count/")
SD_Sham7d <- Read10X(data.dir = "Data/SD-Sham7d_count/")
SD_1h <- Read10X(data.dir = "Data/SD-1h_count/")
SD_4h <- Read10X(data.dir = "Data/SD-4h_count/")
SD_12h <- Read10X(data.dir = "Data/SD-12h_count/")
SD_3d <- Read10X(data.dir = "Data/SD-3d_count/")
SD_7d <- Read10X(data.dir = "Data/SD-7d_count/")

#
SD_Sham <- CreateSeuratObject(counts = SD_Sham, project = "SD_Sham")
SD_Sham3d <- CreateSeuratObject(counts = SD_Sham3d, project = "SD_Sham3d")
SD_Sham7d <- CreateSeuratObject(counts = SD_Sham7d, project = "SD_Sham7d")
SD_1h <- CreateSeuratObject(counts = SD_1h, project = "SD_1h")
SD_4h <- CreateSeuratObject(counts = SD_4h, project = "SD_4h")
SD_12h <- CreateSeuratObject(counts = SD_12h, project = "SD_12h")
SD_3d <- CreateSeuratObject(counts = SD_3d, project = "SD_3d")
SD_7d <- CreateSeuratObject(counts = SD_7d, project = "SD_7d")
####
#SD_Sham
sss=read.table(file = "Data/mart_exportNA.txt", header = T, sep = "\t")
bbb <- row.names(SD_Sham@assays[["RNA"]]@counts)
ccc <- data.frame(feature=setdiff(bbb,sss$Gene.stable.ID))
row.names(ccc) <- ccc$feature
row.names(sss) <- sss$Gene.stable.ID
eee <- data.frame(feature=sss[,2])
row.names(eee) <- sss$Gene.stable.ID
ddd <- rbind(ccc,eee)
loc = match(rownames(SD_Sham@assays$RNA@counts),rownames(ddd))
ddd = ddd[loc,]
row.names(SD_Sham@assays$RNA@counts)=ddd
#SD_Sham3d
sss=read.table(file = "Data/mart_exportNA.txt", header = T, sep = "\t")
bbb <- row.names(SD_Sham3d@assays[["RNA"]]@counts)
ccc <- data.frame(feature=setdiff(bbb,sss$Gene.stable.ID))
row.names(ccc) <- ccc$feature
row.names(sss) <- sss$Gene.stable.ID
eee <- data.frame(feature=sss[,2])
row.names(eee) <- sss$Gene.stable.ID
ddd <- rbind(ccc,eee)
loc = match(rownames(SD_Sham3d@assays$RNA@counts),rownames(ddd))
ddd = ddd[loc,]
row.names(SD_Sham3d@assays$RNA@counts)=ddd
#SD_Sham7d
sss=read.table(file = "Data/mart_exportNA.txt", header = T, sep = "\t")
bbb <- row.names(SD_Sham7d@assays[["RNA"]]@counts)
ccc <- data.frame(feature=setdiff(bbb,sss$Gene.stable.ID))
row.names(ccc) <- ccc$feature
row.names(sss) <- sss$Gene.stable.ID
eee <- data.frame(feature=sss[,2])
row.names(eee) <- sss$Gene.stable.ID
ddd <- rbind(ccc,eee)
loc = match(rownames(SD_Sham7d@assays$RNA@counts),rownames(ddd))
ddd = ddd[loc,]
row.names(SD_Sham7d@assays$RNA@counts)=ddd
#SD_1h
sss=read.table(file = "Data/mart_exportNA.txt", header = T, sep = "\t")
bbb <- row.names(SD_1h@assays[["RNA"]]@counts)
ccc <- data.frame(feature=setdiff(bbb,sss$Gene.stable.ID))
row.names(ccc) <- ccc$feature
row.names(sss) <- sss$Gene.stable.ID
eee <- data.frame(feature=sss[,2])
row.names(eee) <- sss$Gene.stable.ID
ddd <- rbind(ccc,eee)
loc = match(rownames(SD_1h@assays$RNA@counts),rownames(ddd))
ddd = ddd[loc,]
row.names(SD_1h@assays$RNA@counts)=ddd
#SD_4h
sss=read.table(file = "Data/mart_exportNA.txt", header = T, sep = "\t")
bbb <- row.names(SD_4h@assays[["RNA"]]@counts)
ccc <- data.frame(feature=setdiff(bbb,sss$Gene.stable.ID))
row.names(ccc) <- ccc$feature
row.names(sss) <- sss$Gene.stable.ID
eee <- data.frame(feature=sss[,2])
row.names(eee) <- sss$Gene.stable.ID
ddd <- rbind(ccc,eee)
loc = match(rownames(SD_4h@assays$RNA@counts),rownames(ddd))
ddd = ddd[loc,]
row.names(SD_4h@assays$RNA@counts)=ddd
#SD_12h
sss=read.table(file = "Data/mart_exportNA.txt", header = T, sep = "\t")
bbb <- row.names(SD_12h@assays[["RNA"]]@counts)
ccc <- data.frame(feature=setdiff(bbb,sss$Gene.stable.ID))
row.names(ccc) <- ccc$feature
row.names(sss) <- sss$Gene.stable.ID
eee <- data.frame(feature=sss[,2])
row.names(eee) <- sss$Gene.stable.ID
ddd <- rbind(ccc,eee)
loc = match(rownames(SD_12h@assays$RNA@counts),rownames(ddd))
ddd = ddd[loc,]
row.names(SD_12h@assays$RNA@counts)=ddd
#SD_3d
sss=read.table(file = "Data/mart_exportNA.txt", header = T, sep = "\t")
bbb <- row.names(SD_3d@assays[["RNA"]]@counts)
ccc <- data.frame(feature=setdiff(bbb,sss$Gene.stable.ID))
row.names(ccc) <- ccc$feature
row.names(sss) <- sss$Gene.stable.ID
eee <- data.frame(feature=sss[,2])
row.names(eee) <- sss$Gene.stable.ID
ddd <- rbind(ccc,eee)
loc = match(rownames(SD_3d@assays$RNA@counts),rownames(ddd))
ddd = ddd[loc,]
row.names(SD_3d@assays$RNA@counts)=ddd
#SD_7d
sss=read.table(file = "Data/mart_exportNA.txt", header = T, sep = "\t")
bbb <- row.names(SD_7d@assays[["RNA"]]@counts)
ccc <- data.frame(feature=setdiff(bbb,sss$Gene.stable.ID))
row.names(ccc) <- ccc$feature
row.names(sss) <- sss$Gene.stable.ID
eee <- data.frame(feature=sss[,2])
row.names(eee) <- sss$Gene.stable.ID
ddd <- rbind(ccc,eee)
loc = match(rownames(SD_7d@assays$RNA@counts),rownames(ddd))
ddd = ddd[loc,]
row.names(SD_7d@assays$RNA@counts)=ddd




###--creat seurat
SD_Sham <- CreateSeuratObject(counts = SD_Sham@assays$RNA@counts, project = "SD_Sham", min.cells = 3, min.features = 200)
SD_Sham3d <- CreateSeuratObject(counts = SD_Sham3d@assays$RNA@counts, project = "SD_Sham3d", min.cells = 3, min.features = 200)
SD_Sham7d <- CreateSeuratObject(counts = SD_Sham7d@assays$RNA@counts, project = "SD_Sham7d", min.cells = 3, min.features = 200)
SD_1h <- CreateSeuratObject(counts = SD_1h@assays$RNA@counts, project = "SD_1h", min.cells = 3, min.features = 200)
SD_4h <- CreateSeuratObject(counts = SD_4h@assays$RNA@counts, project = "SD_4h", min.cells = 3, min.features = 200)
SD_12h <- CreateSeuratObject(counts = SD_12h@assays$RNA@counts, project = "SD_12h", min.cells = 3, min.features = 200)
SD_3d <- CreateSeuratObject(counts = SD_3d@assays$RNA@counts, project = "SD_3d", min.cells = 3, min.features = 200)
SD_7d <- CreateSeuratObject(counts = SD_7d@assays$RNA@counts, project = "SD_7d", min.cells = 3, min.features = 200)

###--Human-MT、Mouse-mt
SD_Sham[["percent.mt"]] <- PercentageFeatureSet(SD_Sham, pattern = "^Mt-")
SD_Sham3d[["percent.mt"]] <- PercentageFeatureSet(SD_Sham3d, pattern = "^Mt-")
SD_Sham7d[["percent.mt"]] <- PercentageFeatureSet(SD_Sham7d, pattern = "^Mt-")
SD_1h[["percent.mt"]] <- PercentageFeatureSet(SD_1h, pattern = "^Mt-")
SD_4h[["percent.mt"]] <- PercentageFeatureSet(SD_4h, pattern = "^Mt-")
SD_12h[["percent.mt"]] <- PercentageFeatureSet(SD_12h, pattern = "^Mt-")
SD_3d[["percent.mt"]] <- PercentageFeatureSet(SD_3d, pattern = "^Mt-")
SD_7d[["percent.mt"]] <- PercentageFeatureSet(SD_7d, pattern = "^Mt-")

###--HB
HB.genes <- c("Hba-a2","Hba-a3","Hbb-b1","Hbb-bs","Hbg1","Hbe1","Hbe2","Hbz")

HB_m <- match(HB.genes, rownames(SD_Sham@assays[["RNA"]]@counts)) 
HB.genes <- rownames(SD_Sham@assays[["RNA"]]@counts)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
SD_Sham[["percent.HB"]]<-PercentageFeatureSet(SD_Sham, features=HB.genes) 

HB_m <- match(HB.genes, rownames(SD_Sham3d@assays[["RNA"]]@counts)) 
HB.genes <- rownames(SD_Sham3d@assays[["RNA"]]@counts)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
SD_Sham3d[["percent.HB"]]<-PercentageFeatureSet(SD_Sham3d, features=HB.genes) 

HB_m <- match(HB.genes, rownames(SD_Sham7d@assays[["RNA"]]@counts)) 
HB.genes <- rownames(SD_Sham7d@assays[["RNA"]]@counts)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
SD_Sham7d[["percent.HB"]]<-PercentageFeatureSet(SD_Sham7d, features=HB.genes) 

HB_m <- match(HB.genes, rownames(SD_1h@assays[["RNA"]]@counts)) 
HB.genes <- rownames(SD_1h@assays[["RNA"]]@counts)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
SD_1h[["percent.HB"]]<-PercentageFeatureSet(SD_1h, features=HB.genes) 

HB_m <- match(HB.genes, rownames(SD_4h@assays[["RNA"]]@counts)) 
HB.genes <- rownames(SD_4h@assays[["RNA"]]@counts)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
SD_4h[["percent.HB"]]<-PercentageFeatureSet(SD_4h, features=HB.genes) 

HB_m <- match(HB.genes, rownames(SD_12h@assays[["RNA"]]@counts)) 
HB.genes <- rownames(SD_12h@assays[["RNA"]]@counts)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
SD_12h[["percent.HB"]]<-PercentageFeatureSet(SD_12h, features=HB.genes) 

HB_m <- match(HB.genes, rownames(SD_3d@assays[["RNA"]]@counts)) 
HB.genes <- rownames(SD_3d@assays[["RNA"]]@counts)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
SD_3d[["percent.HB"]]<-PercentageFeatureSet(SD_3d, features=HB.genes) 

HB_m <- match(HB.genes, rownames(SD_7d@assays[["RNA"]]@counts)) 
HB.genes <- rownames(SD_7d@assays[["RNA"]]@counts)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
SD_7d[["percent.HB"]]<-PercentageFeatureSet(SD_7d, features=HB.genes) 

###--ribosomes
rownames(SD_Sham)[grepl('^Rp[sl]',rownames(SD_Sham),ignore.case = T)]
rb.genes <- rownames(SD_Sham)[grep("^RP[SL]",rownames(SD_Sham),ignore.case = T)]
C <- GetAssayData(object = SD_Sham, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
SD_Sham <- AddMetaData(SD_Sham, percent.ribo, col.name = "percent.ribo")

rownames(SD_Sham3d)[grepl('^Rp[sl]',rownames(SD_Sham3d),ignore.case = T)]
rb.genes <- rownames(SD_Sham3d)[grep("^RP[SL]",rownames(SD_Sham3d),ignore.case = T)]
C <- GetAssayData(object = SD_Sham3d, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
SD_Sham3d <- AddMetaData(SD_Sham3d, percent.ribo, col.name = "percent.ribo")

rownames(SD_Sham7d)[grepl('^Rp[sl]',rownames(SD_Sham7d),ignore.case = T)]
rb.genes <- rownames(SD_Sham7d)[grep("^RP[SL]",rownames(SD_Sham7d),ignore.case = T)]
C <- GetAssayData(object = SD_Sham7d, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
SD_Sham7d <- AddMetaData(SD_Sham7d, percent.ribo, col.name = "percent.ribo")

rownames(SD_1h)[grepl('^Rp[sl]',rownames(SD_1h),ignore.case = T)]
rb.genes <- rownames(SD_1h)[grep("^RP[SL]",rownames(SD_1h),ignore.case = T)]
C <- GetAssayData(object = SD_1h, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
SD_1h <- AddMetaData(SD_1h, percent.ribo, col.name = "percent.ribo")

rownames(SD_4h)[grepl('^Rp[sl]',rownames(SD_4h),ignore.case = T)]
rb.genes <- rownames(SD_4h)[grep("^RP[SL]",rownames(SD_4h),ignore.case = T)]
C <- GetAssayData(object = SD_4h, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
SD_4h <- AddMetaData(SD_4h, percent.ribo, col.name = "percent.ribo")

rownames(SD_12h)[grepl('^Rp[sl]',rownames(SD_12h),ignore.case = T)]
rb.genes <- rownames(SD_12h)[grep("^RP[SL]",rownames(SD_12h),ignore.case = T)]
C <- GetAssayData(object = SD_12h, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
SD_12h <- AddMetaData(SD_12h, percent.ribo, col.name = "percent.ribo")

rownames(SD_3d)[grepl('^Rp[sl]',rownames(SD_3d),ignore.case = T)]
rb.genes <- rownames(SD_3d)[grep("^RP[SL]",rownames(SD_3d),ignore.case = T)]
C <- GetAssayData(object = SD_3d, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
SD_3d <- AddMetaData(SD_3d, percent.ribo, col.name = "percent.ribo")

rownames(SD_7d)[grepl('^Rp[sl]',rownames(SD_7d),ignore.case = T)]
rb.genes <- rownames(SD_7d)[grep("^RP[SL]",rownames(SD_7d),ignore.case = T)]
C <- GetAssayData(object = SD_7d, slot = "counts")
percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
SD_7d <- AddMetaData(SD_7d, percent.ribo, col.name = "percent.ribo")
###
SD_Sham <- subset(SD_Sham, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
SD_Sham3d <- subset(SD_Sham3d, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
SD_Sham7d <- subset(SD_Sham7d, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
SD_1h <- subset(SD_1h, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
SD_4h <- subset(SD_4h, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
SD_12h <- subset(SD_12h, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
SD_3d <- subset(SD_3d, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
SD_7d <- subset(SD_7d, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)

myfunction1 <- function(xxxx){
  xxxx <- NormalizeData(xxxx, normalization.method = "LogNormalize", scale.factor = 10000)
  xxxx <- FindVariableFeatures(xxxx, selection.method = "vst", nfeatures = 2000)
  return(xxxx)   }

SD_Sham <- myfunction1(SD_Sham)
SD_Sham3d <- myfunction1(SD_Sham3d)
SD_Sham7d <- myfunction1(SD_Sham7d)
SD_1h <- myfunction1(SD_1h)
SD_4h <- myfunction1(SD_4h)
SD_12h <- myfunction1(SD_12h)
SD_3d <- myfunction1(SD_3d)
SD_7d <- myfunction1(SD_7d)

AB.anchors <- FindIntegrationAnchors(object.list = list(SD_Sham,SD_Sham3d,SD_Sham7d,SD_1h,SD_4h,SD_12h,SD_3d,SD_7d), anchor.features = 2000,dims = 1:30)  
AB <- IntegrateData(anchorset = AB.anchors, dims = 1:30)

DefaultAssay(AB) <- "integrated"
AB <- ScaleData(AB, features = rownames(AB))
AB <- RunPCA(AB, features = VariableFeatures(object = AB))
print(AB[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(AB, dims = 1:2, reduction = "pca")
DimPlot(AB, reduction = "pca")
DimHeatmap(AB, dims = 1:15, cells = 500, balanced = TRUE)

AB <- JackStraw(AB, num.replicate = 100)
AB <- ScoreJackStraw(AB, dims = 1:20)

plot1 <- JackStrawPlot(AB, dims = 1:20)
ggsave(filename = "JackStrawPlot.pdf", plot = plot1, device = 'pdf', width = 14, height = 12, units = 'cm')
plot2 <- ElbowPlot(AB)
ggsave(filename = "ElbowPlot.pdf", plot = plot2, device = 'pdf', width = 14, height = 12, units = 'cm')
rm('plot1','plot2')
save(AB,file="AB.RData")

AB <- FindNeighbors(AB, dims = 1:30)
AB <- FindClusters(AB, resolution = 1.2)
AB <- RunUMAP(AB, dims = 1:30)
DimPlot(AB, reduction = "umap")
#AB <- RunTSNE(AB, dims = 1:30)
#DimPlot(AB, reduction = "tsne")

#
table(AB@meta.data$orig.ident)
table(AB@meta.data$seurat_clusters)
table(AB@meta.data$celltype)
AB$age=str_replace(AB$orig.ident,"SD_","")#?????

AB@meta.data[["age"]]<-factor(AB@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
DefaultAssay(AB) <- "RNA"
save(AB,file="AB_new.RData")

p7 <- DimPlot(AB, reduction = "umap", group.by = "age", pt.size=0.01,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2_Cluster.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

DefaultAssay(AB) <- "RNA"
#
markers <- FindAllMarkers(AB, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox") 
write.table(markers,file="markers.txt",quote=F,sep="\t",row.names=F,col.names=T)

AB100=subset(AB,downsample=100)
markers100 <- FindAllMarkers(AB100, logfc.threshold = 0.25, min.pct = 0.5, only.pos = F, test.use = "wilcox") 
write.table(markers100,file="markers100.txt",quote=F,sep="\t",row.names=F,col.names=T)

#
plotPig_0 <- VlnPlot(AB, group.by = "age",features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB", "percent.ribo"), cols = c("#bE555D","#FFE76F","#72535F","#98D6E6","#80C684",'#DB6626'), pt.size = 0.001 ,ncol = 5)
ggsave("ALL_control.pdf", width = 78, height = 25, units = "cm")

jjplot <- jjDotPlot(object = AB,
                    gene = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1","Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Nrp2","Ntf3","Ntng1","Sorcs3"),
                    gene.order = c("Hba-a2","Hbb-b1","Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Cspg4","Pdgfra","Inpp5d","Csf1r","Lgmn","P2ry12","Ifngr1","Bcas1","Nfasc","Sema5a","Tcf7l2","Aqp4","Gfap","Slc1a3","Slc1a2","Myt1","Col1a1","Col1a2","Igfbp7","Pdgfrb","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Sox6","Olig1","Olig2","Mbp","Mobp","Plp1","Egfr","Ascl1","Cd163","Bdnf","Fos","Jun","Nr4a3","Met","Ube2h","Mki67","Top2a","Nrp2","Ntf3","Ntng1","Sorcs3"),
                    #cluster.order = c(7,4,0,6,2,3,1,5,8,10,9),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 7)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot.pdf", plot = jjplot, device = 'pdf', width = 40, height = 40, units = 'cm')


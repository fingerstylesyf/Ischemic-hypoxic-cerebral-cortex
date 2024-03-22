#-----------------Astrocytes-----------------#
library(multtest)
library(Seurat)
library(dplyr)
library(mindr)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)
library(scRNAtoolVis)
#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
library(Scillus)
#devtools::install_github('junjunlab/scRNAtoolVis')
#remotes::install_github("lyc-1995/MySeuratWrappers")
#devtools::install_github("xmc811/Scillus", ref = "development")

#-----------data integrating------------#
aging.list<-SplitObject(AB, split.by = "celltype")
scRNA <- aging.list$AST
save(scRNA,file="AST.RData")
#
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
#
scRNA <- FindNeighbors(scRNA, dims = 1:12)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:12)
DimPlot(scRNA, reduction = "umap")
#scRNA <- RunTSNE(scRNA, dims = 1:20)
#DimPlot(scRNA, reduction = "tsne")
#
p7 <- DimPlot(scRNA, reduction = "umap", group.by = "age", pt.size=0.1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", pt.size=0.7, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap1 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap1.pdf", plot = umap1, device = 'pdf', width = 34, height = 14.5, units = 'cm')
# Delete unwanted cells
scRNA <- scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("0",	"1",	"2",	"3",	"4",	"5", "6",	"7",	"8",	"9",	"10",	"11",	"12",	"13")]
#
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
#
scRNA <- FindNeighbors(scRNA, dims = 1:11)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims = 1:11)
DimPlot(scRNA, reduction = "umap")
#
p7 <- DimPlot(scRNA, reduction = "umap", group.by = "age", pt.size=0.1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", pt.size=0.7, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2.pdf", plot = umap2, device = 'pdf', width = 34, height = 14.5, units = 'cm')
#
scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$seurat_clusters)
table(scRNA@meta.data$celltype)
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["age"]])
DefaultAssay(scRNA) <- "RNA"
save(scRNA,file="scRNA.RData")
###--
#
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.5, only.pos = F, test.use = "wilcox") 
write.table(markers,file="markers-cluster.txt",quote=F,sep="\t",row.names=F,col.names=T)
#
jjplot <- jjDotPlot(object = scRNA,
                    gene = c("Slc1a2","Slc1a3","Aqp4","Gfap","Angpt1","Egfr","S100a10","Cd44","Plaur","S1pr3","Emp1","Disp3","Alk","Sema5a","Jun","Fos","Dnah9","Dnah3","Top2a","Cenpf","Birc5","Rrm2","Hells","Rbfox3","Dcx","Ptx3","Tfap2c","Hopx","Sox2","Ascl1","Egr1","Dlx1","Gad1","Bcl11a","Fabp7","Hes5","Nr2e1","Aldh1l1","Gjb6","Cdk2",'Cd9',"Ifngr1","Id3","Id2","Sox9","Lfng","Notch1","Notch2","Bmp6","Tmem176a","Prom1","Cd24","Tgfbr1","Sema5b","S100b","Stat3","Igf1r","Tnf","Ngf","Ltbp1","Sbno2","Nfatc2","Dclk1","Bcl2","Hspa4","Nfkb1","Jak2","Mapk1","Olig2","Socs3","Prkaca","Cntf","Mapk14","Ccl2","Ccl3","Ccl4","Il6","Bdnf","Gdnf","Fgf2","Vegfa"),
                    gene.order = c("Slc1a2","Slc1a3","Aqp4","Gfap","Angpt1","Egfr","S100a10","Cd44","Plaur","S1pr3","Emp1","Disp3","Alk","Sema5a","Jun","Fos","Dnah9","Dnah3","Top2a","Cenpf","Birc5","Rrm2","Hells","Rbfox3","Dcx","Ptx3","Tfap2c","Hopx","Sox2","Ascl1","Egr1","Dlx1","Gad1","Bcl11a","Fabp7","Hes5","Nr2e1","Aldh1l1","Gjb6","Cdk2",'Cd9',"Ifngr1","Id3","Id2","Sox9","Lfng","Notch1","Notch2","Bmp6","Tmem176a","Prom1","Cd24","Tgfbr1","Sema5b","S100b","Stat3","Igf1r","Tnf","Ngf","Ltbp1","Sbno2","Nfatc2","Dclk1","Bcl2","Hspa4","Nfkb1","Jak2","Mapk1","Olig2","Socs3","Prkaca","Cntf","Mapk14","Ccl2","Ccl3","Ccl4","Il6","Bdnf","Gdnf","Fgf2","Vegfa"),
                    cluster.order = 0:8,
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)
jjplot
ggsave(filename = "jjDotPlot.pdf", plot = jjplot, device = 'pdf', width = 28.5, height = 10.5, units = 'cm')
###--rename cluster
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8")
new.cluster.ids <- c("AST1","AST2","AST3","AST5","AST4","AST7","AST8","AST6","EPC")
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$celltype)
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["age"]])
#
Idents(scRNA) <- factor(scRNA$celltype, levels = c("AST1","AST2","AST3","AST4","AST5","AST6","AST7","AST8","EPC"))
DimPlot(scRNA, reduction = "umap", label = TRUE)
save(scRNA,file="AST_celltype.RData")
#
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox")  ##耗时久
write.table(markers,file="markers-celltype.txt",quote=F,sep="\t",row.names=F,col.names=T)

################################################################
###-----------plot------------####
#---------UMAP(celltype)
p7 <- DimPlot(AB, reduction = "umap", group.by = "age", pt.size=0.01,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(AB, reduction = "umap", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE,raster=FALSE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap_celltype.pdf", plot = umap2, device = 'pdf', width = 43, height = 17, units = 'cm')

#---------Dotplot(celltype)
jjplot <- jjDotPlot(object = scRNA,
                    gene = c("Slc1a3","Gfap","Cd44","Emp1","Sbno2","Dclk1","Hmga1","Hmga2","Camk2d","Fosl1","Nrg1","Jun","Fos","Top2a","Cenpf","Egfr","Dnah9","Dnah3"),
                    gene.order = c("Slc1a3","Gfap","Cd44","Dclk1","Camk2d","Emp1","Sbno2","Hmga1","Hmga2","Fosl1","Nrg1","Jun","Fos","Top2a","Cenpf","Egfr","Dnah9","Dnah3"),
                    id = 'celltype',
                    cluster.order = c("AST1",	"AST2",	"AST3",	"AST4",	"AST5",	"AST6",	"AST7",	"AST8",	"EPC"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot_celltype.pdf", plot = jjplot, device = 'pdf', width = 25, height = 11, units = 'cm')

#-----------Mfuzz----------#
library(Mfuzz)
library(Seurat)
library(progeny)
library(tidyr)
library(tibble)
library(dplyr)
library(viridisLite)
library(ggplot2)
#BiocManager::install("progeny")
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
#
scRNA=subset(scRNA,features= rownames(scRNA@assays$RNA@scale.data))
scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
scRNA=ScaleData(scRNA)
##
age.averages <- AverageExpression(scRNA,group.by = "age")
head(age.averages[["RNA"]])
head(age.averages[["integrated"]])
df1 <- age.averages[["RNA"]]###这里的RNA就是intergrated
#
mat <- as.matrix(df1)
dt <- new("ExpressionSet",exprs = mat)
dim(dt)
#
dt.r <- filter.NA(dt, thres=0.25)
dt.f <- fill.NA(dt.r,mode="mean")
tmp <- filter.std(dt.f,min.std=0.15)
dt.s <- standardise(tmp)
df.s <- dt.s@assayData$exprs
head(df.s)
#
m1 <- mestimate(dt.s)
#
set.seed(007)
cl <- mfuzz(dt.s,c=12,m=m1)
#
mfuzz.plot(dt.s,cl,mfrow=c(4,4), new.window= FALSE, time.labels=colnames(dt.s))
#
library(RColorBrewer)
mycol <- c("cyan","yellow","orangered")
mycol <- c("#737fb4","#168c90","#f7c99b","#421c77")
mycolor <- colorRampPalette(mycol)(20)
#
pdf("AST expression changes trends_1.pdf", width = 14,height = 13)
mfuzz.plot(dt.s,cl,mfrow=c(4,4), new.window= FALSE, time.labels=colnames(dt.s), colo = mycolor)
dev.off()
#
gene_cluster <- data.frame(df.s,cluster=cl$cluster)
#
write.csv(gene_cluster, 'Mfuzz genes set.csv')

#-----------Mfuzz genes set score----------#
###
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix)
library(readxl)
library(ggpubr)
#
HIE_gene<-read_xlsx("Mfuzz genes set.xlsx")
gene<-as.list(HIE_gene)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "HIE")
colnames(AB@meta.data)[10]<-"HIE_Score"
#####
AB@meta.data[["celltype"]]<-factor(AB@meta.data[["celltype"]], levels=c("AST1","AST2","AST3","AST4","AST5","AST6","AST7","AST8","EPC"))


my_comparisons <- list(c('0','5'),c('1','2'))
my9color <- c("#fbd15c","#6F7796","#AB4361","#71CB89","#Afd2ce","#2b5aab","#03776a","#c7ad6c","#0EB8C9")
library(ggpubr)
p.AddModuleScore <- ggviolin(AB@meta.data, x = "celltype", y = "HIE_Score",
                             color = "celltype",add = 'mean_sd',fill = 'celltype',
                             add.params = list(color = "black")) + 
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color) +
  #theme(axis.text.x.bottom = element_text(angle = 0,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')
p.AddModuleScore
ggsave('Mfuzz genes(AST4)_vlo.pdf',width = 5,height = 4)


HIE_gene<-read_xlsx("Mfuzz genes set.xlsx")
gene<-as.list(HIE_gene)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "HIE")
colnames(AB@meta.data)[10]<-"HIE_Score"
#####
AB@meta.data[["celltype"]]<-factor(AB@meta.data[["celltype"]], levels=c("AST1","AST2","AST3","AST4","AST5","AST6","AST7","AST8","EPC"))


my_comparisons <- list(c('0','5'),c('1','2'))
my9color <- c("#fbd15c","#6F7796","#AB4361","#71CB89","#Afd2ce","#2b5aab","#03776a","#c7ad6c","#0EB8C9")
library(ggpubr)
p.AddModuleScore <- ggviolin(AB@meta.data, x = "celltype", y = "HIE_Score",
                             color = "celltype",add = 'mean_sd',fill = 'celltype',
                             add.params = list(color = "black")) + 
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color) +
  #theme(axis.text.x.bottom = element_text(angle = 0,vjust = 0.5,hjust = 1)) + 
  NoLegend() + labs(x = '')
p.AddModuleScore
ggsave('Mfuzz genes(AST6)_vlo.pdf',width = 5,height = 4)


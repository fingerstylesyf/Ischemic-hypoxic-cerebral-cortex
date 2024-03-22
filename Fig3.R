#-----------------Microglia-----------------#
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
scRNA <- aging.list$MG
save(scRNA,file="MG.RData")
#
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
plot2 <- ElbowPlot(scRNA)
plot2
#
scRNA <- FindNeighbors(scRNA, dims = 1:15)
scRNA <- FindClusters(scRNA, resolution = 0.5)
scRNA <- RunUMAP(scRNA, dims = 1:15)
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
scRNA <- scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("0",	"1",	"2",	"3",	"4",	"5", "6",	"7",	"8",	"10")]
#
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
plot2 <- ElbowPlot(scRNA)
plot2
save(scRNA,file="scRNA.RData")
#
scRNA <- FindNeighbors(scRNA, dims = 1:12)
scRNA <- FindClusters(scRNA, resolution = 0.5)
scRNA <- RunUMAP(scRNA, dims = 1:12)
DimPlot(scRNA, reduction = "umap")
#
p7 <- DimPlot(scRNA, reduction = "umap", group.by = "age", pt.size=0.1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", pt.size=0.7, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2.pdf", plot = umap2, device = 'pdf', width = 34, height = 14.5, units = 'cm')
#
#scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$seurat_clusters)
table(scRNA@meta.data$celltype)
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["age"]])
DefaultAssay(scRNA) <- "RNA"
save(scRNA,file="scRNA2.RData")
###--
#
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox") 
write.table(markers,file="markers-cluster.txt",quote=F,sep="\t",row.names=F,col.names=T)
#
jjplot <- jjDotPlot(object = scRNA,
                    gene = c("Csf1r","P2ry12","Cx3cr1","Dab2","Mrc1","Cd163","Apoe","Trem2","Cd68","Tyrobp","Lrp4","Ldlrad4","Abca1","Apobec1","Axl","Cd47","Itgam","Ccr5","Il18","Spp1","Arg1","Piezo1","Igf1","Smad3","Top2a","Cenpf","Rbfox3","Satb2","Gad1","Gfap","Slc1a3","Pdgfra","Cspg4","Entpd1","Plau","S100b","Gpnmb"),
                    gene.order = c("Csf1r","P2ry12","Cx3cr1","Dab2","Mrc1","Cd163","Apoe","Trem2","Cd68","Tyrobp","Lrp4","Ldlrad4","Abca1","Apobec1","Axl","Cd47","Itgam","Ccr5","Il18","Spp1","Arg1","Piezo1","Igf1","Smad3","Top2a","Cenpf","Rbfox3","Satb2","Gad1","Gfap","Slc1a3","Pdgfra","Cspg4","Entpd1","Plau","S100b","Gpnmb"),
                    cluster.order = 0:9,
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)
jjplot
ggsave(filename = "jjDotPlot.pdf", plot = jjplot, device = 'pdf', width = 28.5, height = 10.5, units = 'cm')
###--rename cluster
current.cluster.ids <- c("0","1","2","3","4","5","6","7","8")
new.cluster.ids <- c("MG1","MG2","MG4","MG8","MG3","MG3","MG5","MG6","MG7")
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$celltype)
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["age"]])
Idents(scRNA) <- factor(scRNA$celltype, levels = c("MG1","MG2","MG3","MG4","MG5","MG6","MG7","MG8"))

DimPlot(scRNA, reduction = "umap", label = TRUE)
save(scRNA,file="MG-celltype.RData")
#
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox")  ##耗时久
write.table(markers,file="markers-celltype.txt",quote=F,sep="\t",row.names=F,col.names=T)


################################################################
###-----------plot------------####
#---------UMAP(celltype)
p7 <- DimPlot(scRNA, reduction = "umap", group.by = "age", pt.size=0.1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.7, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2-celltype.pdf", plot = umap2, device = 'pdf', width = 34, height = 14.5, units = 'cm')

#---------Dotplot(celltype)
jjplot <- jjDotPlot(object = scRNA,
                    gene = c("Csf1r","P2ry12","Itgam","Ccr5","Il18","Apoe","Abca1","Trem2","S100b","Dab2","Mrc1","Cd163","Itgax","Spp1","Gpnmb","Piezo1","Igf1","Csf1","Plau","Olr1","Top2a","Cenpf","Tmem119","Lgals3","Cd63","Cd9","Nes","Vat1","Ccl3","Pik3ip1","Fosl1","Bach1"),
                    gene.order = c("Csf1r","P2ry12","Itgam","Ccr5","Il18","Apoe","Abca1","Trem2","S100b","Dab2","Mrc1","Cd163","Spp1","Gpnmb","Piezo1","Igf1","Csf1","Plau","Olr1","Itgax","Top2a","Cenpf","Tmem119","Lgals3","Cd63","Cd9","Nes","Vat1","Ccl3","Pik3ip1","Fosl1","Bach1"),
                    id = 'celltype',
                    cluster.order = c("MG1",	"MG2",	"MG3",	"MG4",	"MG5",	"MG6",	"MG7",	"MG8"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "jjDotPlot-celltype.pdf", plot = jjplot, device = 'pdf', width = 24.5, height = 11, units = 'cm')

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
##
#
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
#
scRNA=subset(scRNA,features= rownames(scRNA@assays$RNA@scale.data))
scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
####
scRNA=ScaleData(scRNA)
age.averages <- AverageExpression(scRNA,group.by = "age")
head(age.averages[["RNA"]])
head(age.averages[["integrated"]])
df1 <- age.averages[["integrated"]]
mat <- as.matrix(df1)
#
dt <- new("ExpressionSet",exprs = mat)
dim(dt)
dt.r <- filter.NA(dt, thres=0.25)
#
dt.f <- fill.NA(dt.r,mode="mean")
tmp <- filter.std(dt.f,min.std=0)
#
dt.s <- standardise(tmp)
df.s <- dt.s@assayData$exprs
head(df.s)
m1 <- mestimate(dt.s)
set.seed(007)
cl <- mfuzz(dt.s,c=6,m=m1)
mfuzz.plot(dt.s,cl,mfrow=c(2,3), new.window= FALSE, time.labels=colnames(dt.s))

library(RColorBrewer)
mycol <- c("cyan","yellow","orangered")
mycol <- c("#737fb4","#168c90","#f7c99b","#421c77")
mycolor <- colorRampPalette(mycol)(20)

#
pdf("MG expression changes trends_1.pdf", width = 14,height = 7)
mfuzz.plot(dt.s,cl,mfrow=c(2,3), new.window= FALSE, time.labels=colnames(dt.s), colo = mycolor)
dev.off()
#
cl$size
gene_cluster <- data.frame(df.s,cluster=cl$cluster)
#
write.csv(gene_cluster, 'Mfuzz genes set.csv')

#-----------Mfuzz genes set score----------#
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix)
library(readxl)
library(ggpubr)

#
HIE_gene<-read_xlsx("Mfuzz genes set score.xlsx")
gene<-as.list(HIE_gene)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "HIE")
colnames(AB@meta.data)[9]<-"HIE_Score"
#####
AB@meta.data[["celltype"]]<-factor(AB@meta.data[["celltype"]], levels=c("MG1","MG2","MG3","MG4","MG5","MG6","MG7","MG8"))


my_comparisons <- list(c('0','5'),c('1','2'))
my9color <- c("#f26f66","#4962a1","#8572ad",'#3082b4',"#539075",'#63b9bd',"#ccaf74","#2F4858")
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
ggsave('Mfuzz genes(MG6)_vlo.pdf',width = 5,height = 4)

HIE_gene<-read_xlsx("Mfuzz genes set score.xlsx")
gene<-as.list(HIE_gene)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "HIE")
colnames(AB@meta.data)[9]<-"HIE_Score"
#####
AB@meta.data[["celltype"]]<-factor(AB@meta.data[["celltype"]], levels=c("MG1","MG2","MG3","MG4","MG5","MG6","MG7","MG8"))

my_comparisons <- list(c('0','5'),c('1','2'))
my9color <- c("#f26f66","#4962a1","#8572ad",'#3082b4',"#539075",'#63b9bd',"#ccaf74","#2F4858")
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
ggsave('Mfuzz genes(MG7)_vlo.pdf',width = 5,height = 4)

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
scRNA <- aging.list$Neuron
save(scRNA,file="Neuron.RData")
#
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
plot2 <- ElbowPlot(scRNA)
plot2
#
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA, resolution = 1)
scRNA <- RunUMAP(scRNA, dims = 1:20)#
DimPlot(scRNA, reduction = "umap")
#scRNA <- RunTSNE(scRNA, dims = 1:20)
#DimPlot(scRNA, reduction = "tsne")
#
p7 <- DimPlot(scRNA, reduction = "umap", group.by = "age", pt.size=0.1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap1 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap1.pdf", plot = umap1, device = 'pdf', width = 34, height = 14.5, units = 'cm')

DefaultAssay(scRNA) <- "RNA"
#
AB300 <- subset(scRNA,downsample=300)
markers <- FindAllMarkers(AB300, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox")
write.table(markers,file="markers-cluster.txt",quote=F,sep="\t",row.names=F,col.names=T)

jjplot <- jjDotPlot(object = scRNA,
                    gene = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Adarb2","Lhx6","Inpp5d","Csf1r","Lgmn","Ptprc","P2ry12","Ifngr1","Bcas1","Sema5a","Tcf7l2","Plp1","Olig1","Olig2","Sox6","Mbp","Mobp","Pdgfra","Myt1","Aqp4","Gfap","Slc1a3","Slc1a2","Col1a1","Col1a2","Igfbp7","Pdgfrb","Rgs5","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Top2a","Cenpf","Kif4a","Fosb","Fos"),
                    gene.order = c("Rbfox3","Snap25","Camk2a","Satb2","Slc17a7","Gad1","Gad2","Adarb2","Lhx6","Inpp5d","Csf1r","Lgmn","Ptprc","P2ry12","Ifngr1","Bcas1","Sema5a","Tcf7l2","Plp1","Olig1","Olig2","Sox6","Mbp","Mobp","Pdgfra","Myt1","Aqp4","Gfap","Slc1a3","Slc1a2","Col1a1","Col1a2","Igfbp7","Pdgfrb","Rgs5","Adgrl4","Vwf","Reln","Dnah9","Dnah3","Top2a","Cenpf","Kif4a","Fosb","Fos"),
                    cluster.order = 0:29,
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)
jjplot
ggsave(filename = "jjDotPlot.pdf", plot = jjplot, device = 'pdf', width = 40.5, height = 22, units = 'cm')

# Delete unwanted cells
scRNA <- scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20",	"21")]
Dying <- scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("16")]
save(Dying,file="Dying.RData")
#
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
plot2 <- ElbowPlot(scRNA,ndims = 30)
plot2
save(scRNA,file="scRNA.RData")
#
scRNA <- FindNeighbors(scRNA, dims = 1:24)
scRNA <- FindClusters(scRNA, resolution = 0.7)
scRNA <- RunUMAP(scRNA, dims = 1:24)
DimPlot(scRNA, reduction = "umap")
#
p7 <- DimPlot(scRNA, reduction = "umap", group.by = "age", pt.size=0.1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2.pdf", plot = umap2, device = 'pdf', width = 34, height = 14.5, units = 'cm')
##
#scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$seurat_clusters)
table(scRNA@meta.data$celltype)
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["age"]])
DefaultAssay(scRNA) <- "RNA"
save(scRNA,file="scRNA2.RData")
#
AB100 <- subset(scRNA,downsample=100)
markers <- FindAllMarkers(AB100, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox")
write.table(markers,file="markers-cluster2.txt",quote=F,sep="\t",row.names=F,col.names=T)
#
jjplot <- jjDotPlot(object = scRNA,
                    gene = c("Rbfox3","Slc17a7","Gad2","Cux2","Rfx3","Rorb","Sulf2","Nrp2","Ntf3","Ntng1","Sema3e","Gnb4","Nr4a2","Bcl11b","Fezf2","Tshz2","Bcl6","Tle4","Foxp2","Otof","Ccn2","Mme","Sst","Prox1","Vip","Sv2c",'Meis2','Pbx3','Tshz1',"Zfhx4","Reln","Ier2","Hspa5","Sod2","Atf4","Fos","Bdnf","Ndnf","Dcx"),
                    gene.order = c("Rbfox3","Slc17a7","Gad2","Cux2","Rfx3","Rorb","Sulf2","Nrp2","Ntf3","Ntng1","Sema3e","Gnb4","Nr4a2","Bcl11b","Fezf2","Tshz2","Bcl6","Tle4","Foxp2","Otof","Ccn2","Mme","Sst","Prox1","Vip","Sv2c",'Meis2','Pbx3','Tshz1',"Zfhx4","Reln","Ier2","Hspa5","Sod2","Atf4","Fos","Bdnf","Ndnf","Dcx"),
                    cluster.order = 0:23,
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)
jjplot
ggsave(filename = "jjDotPlot2.pdf", plot = jjplot, device = 'pdf', width = 40.5, height = 22, units = 'cm')

# Delete unwanted cells
scRNA <- scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"17",	"18",	"19",	"20",	"22",	"23")]
#
DefaultAssay(scRNA) <- "integrated"
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
plot2 <- ElbowPlot(scRNA,ndims = 30)
plot2
save(scRNA,file="scRNA3.RData")
#
scRNA <- FindNeighbors(scRNA, dims = 1:22)
scRNA <- FindClusters(scRNA, resolution = 0.6)
scRNA <- RunUMAP(scRNA, dims = 1:22)
DimPlot(scRNA, reduction = "umap")

#
p7 <- DimPlot(scRNA, reduction = "umap", group.by = "age", pt.size=0.1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "ident", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap3.pdf", plot = umap2, device = 'pdf', width = 34, height = 14.5, units = 'cm')

##
#scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
table(scRNA@meta.data$orig.ident)
table(scRNA@meta.data$seurat_clusters)
table(scRNA@meta.data$celltype)
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["age"]])
DefaultAssay(scRNA) <- "RNA"
save(scRNA,file="scRNA4.RData")
###--
#
AB100 <- subset(scRNA,downsample=100)
markers <- FindAllMarkers(AB100, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox") 
write.table(markers,file="markers-cluster3.txt",quote=F,sep="\t",row.names=F,col.names=T)
#
jjplot <- jjDotPlot(object = scRNA,
                    gene = c("Rbfox3","Slc17a7","Gad2","Cux2","Rfx3","Rorb","Sulf2","Nrp2","Ntf3","Ntng1","Sema3e","Gnb4","Nr4a2","Bcl11b","Fezf2","Tshz2","Bcl6","Tle4","Foxp2","Otof","Ccn2","Mme","Sst","Prox1","Vip","Sv2c",'Meis2','Pbx3','Tshz1',"Zfhx4","Reln","Ier2","Hspa5","Sod2","Atf4","Fos","Bdnf","Ndnf","Dcx","Nrp2"),
                    gene.order = c("Rbfox3","Slc17a7","Gad2","Cux2","Rfx3","Rorb","Sulf2","Nrp2","Ntf3","Ntng1","Sema3e","Gnb4","Nr4a2","Bcl11b","Fezf2","Tshz2","Bcl6","Tle4","Foxp2","Otof","Ccn2","Mme","Sst","Prox1","Vip","Sv2c",'Meis2','Pbx3','Tshz1',"Zfhx4","Reln","Ier2","Hspa5","Sod2","Atf4","Fos","Bdnf","Ndnf","Dcx","Nrp2"),
                    cluster.order = 0:20,
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)
jjplot
ggsave(filename = "jjDotPlot3.pdf", plot = jjplot, device = 'pdf', width = 40.5, height = 22, units = 'cm')
#
current.cluster.ids <- c("0",	"1",	"2",	"3",	"4",	"5",	"6",	"7",	"8",	"9",	"10",	"11",	"12",	"13",	"14",	"15",	"16",	"17",	"18",	"19",	"20")
new.cluster.ids <- c("L6 CT",	"L2-3 IT",	"L2-3 IT",	"L6 IT",	"L4 IT",	"L2-3 IT",	"Mme",	"L4 IT",	"L2-3 IT",	"Sst",	"Vip",	"L5 ET",	"L5 IT",	"Sv2c",	"Meis2",	"L2-3 IT",	"L5 NP",	"L6 CLA",	"L6b CT",	"L6 CT",	"CR")
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$celltype)
table(scRNA@meta.data[["celltype"]], scRNA@meta.data[["age"]])
Idents(scRNA) <- factor(scRNA$celltype, levels = c("CR",	"L2-3 IT",	"L4 IT",	"L5 IT",	"L6 IT",	"L5 NP",	"L5 ET",	"L6 CLA",	"L6 CT",	"L6b CT",	"Meis2",	"Vip",	"Sv2c",	"Sst",	"Mme"))

DimPlot(scRNA, reduction = "umap", label = TRUE)
save(scRNA,file="scRNA-celltype.RData")

AB100 <- subset(scRNA,downsample=100)
markers <- FindAllMarkers(AB100, logfc.threshold = 0.25, min.pct = 0.25, only.pos = F, test.use = "wilcox")
write.table(markers,file="markers-celltype.txt",quote=F,sep="\t",row.names=F,col.names=T)

################################################################
###-----------plot------------####
#---------UMAP(celltype)
p7 <- DimPlot(scRNA, reduction = "umap", group.by = "age", pt.size=0.1)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap3.pdf", plot = umap2, device = 'pdf', width = 34, height = 14.5, units = 'cm')

#---------ExN_Dotplot
ExN <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("L6b CT",	"L6 CT",	"L6 CLA",	"L5 ET",	"L5 NP",	"L6 IT",	"L5 IT",	"L4 IT",	"L2-3 IT")]

jjplot <- jjDotPlot(object = ExN,
                    gene = c("Rbfox3","Slc17a7","Satb2","Cux2","Rorb","Sulf2","Sema3e","Gnb4","Nr4a2","Bcl11b","Fezf2","Tshz2","Bcl6","Tle4","Foxp2","Otof","Ccn2"),
                    gene.order = c("Rbfox3","Satb2","Slc17a7","Cux2","Rorb","Sulf2","Sema3e","Gnb4","Nr4a2","Bcl11b","Fezf2","Tshz2","Bcl6","Tle4","Foxp2","Otof","Ccn2"),
                    id = 'celltype',
                    cluster.order = c("L6b CT",	"L6 CT",	"L5 ET",	"L5 NP", "L6 CLA",	"L6 IT",	"L5 IT",	"L4 IT",	"L2-3 IT"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "ExN-Dotplot.pdf", plot = jjplot, device = 'pdf', width = 21.5, height = 15, units = 'cm')

#---------InN_Dotplot
jjplot <- jjDotPlot(object = InN,
                    gene = c("Nxph1",	"Erbb4",	"Kcnc2",	"Galntl6",	"Sox6",	"Adamts8",	"Adamts15",	"Maf",	"Kcnip1",	"Igsf11",	"Zfp536",	"Ptchd4",	"Rps6ka5",	"Eya4",	"Cntnap2",	"Ptprm",	"Ppargc1a",	"Cdh9",	"Slit2",	"Gria1",	"6330411D24Rik",	"Gabrg3",	"Kazn",	"Atrnl1",	"Mme",	"Dlgap2",	"Thsd7a",	"Fam13c",	"Sh3rf3",	"Pde5a",	"Mkx",	"Etl4",	"Zfp804b",	"Gad2",	"Il1rapl1",	"Zfp804a",	"Tafa2",	"Spock3",	"Btbd11",	"Ankrd55",	"Klf12",	"Epha6",	"Mgat4c",	"Zfp385d",	"Dscaml1",	"Tmem132c",	"Specc1",	"Limch1",	"Adamts17",	"Sntg1",	"Grip1",	"Ptprt",	"Rora",	"Gria4",	"Npas3",	"Cntnap5c",	"Runx2",	"Phf20l1",	"Ubash3b",	"Luzp2",	"Dgkg",	"AABR07007642.1",	"Elavl2",	"Cntn6",	"Samd5",	"Grm3",	"Mef2a",	"St6gal2",	"St8sia4",	"Cacna2d2",	"Myo1b",	"Dlgap1",	"Zswim6",	"Grm5",	"Fat3",	"Kcns3",	"Maml3",	"Slc44a5",	"Sash1",	"Vwc2"),
                    gene.order = c("Nxph1",	"Erbb4",	"Kcnc2",	"Galntl6",	"Sox6",	"Grik1",	"Rbms3",	"Maf",	"Kcnip1",	"Igsf11",	"Zfp536",	"Ptchd4",	"Rps6ka5",	"Eya4",	"Cntnap2",	"Ptprm",	"Ppargc1a",	"Cdh9",	"Slit2",	"Gria1",	"6330411D24Rik",	"Gabrg3",	"Kazn",	"Atrnl1",	"Mme",	"Dlgap2",	"Thsd7a",	"Fam13c",	"Sh3rf3",	"Pde5a",	"Mkx",	"Etl4",	"Zfp804b",	"Gad2",	"Il1rapl1",	"Zfp804a",	"Tafa2",	"Spock3",	"Btbd11",	"Ankrd55",	"Klf12",	"Epha6",	"Mgat4c",	"Zfp385d",	"Dscaml1",	"Tmem132c",	"Specc1",	"Limch1",	"Adamts17",	"Sntg1",	"Grip1",	"Ptprt",	"Rora",	"Gria4",	"Npas3",	"Cntnap5c",	"Runx2",	"Phf20l1",	"Ubash3b",	"Luzp2",	"Dgkg",	"AABR07007642.1",	"Elavl2",	"Cntn6",	"Samd5",	"Grm3",	"Mef2a",	"St6gal2",	"St8sia4",	"Cacna2d2",	"Myo1b",	"Dlgap1",	"Zswim6",	"Grm5",	"Fat3",	"Kcns3",	"Maml3",	"Slc44a5",	"Sash1",	"Vwc2"),
                    id = 'celltype',
                    cluster.order = c("CR",	"Meis2",	"Sv2c",	"Vip",	"Sst",	"Mme"),
                    ytree = F,
                    rescale = T,
                    rescale.min = 0,
                    rescale.max = 2,dot.max = 8)#+ggplot2:::coord_flip()
jjplot
ggsave(filename = "InN-Dotplot.pdf", plot = jjplot, device = 'pdf', width = 19.5, height = 12.5, units = 'cm')

#-----------Mfuzz(IT)----------#
library(Mfuzz)
library(Seurat)
library(progeny)
library(tidyr)
library(tibble)
library(dplyr)
library(viridisLite)
library(ggplot2)
#
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
#
scRNA=subset(scRNA,features= rownames(scRNA@assays$RNA@scale.data))
scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
scRNA=ScaleData(scRNA)
#
age.averages <- AverageExpression(scRNA,group.by = "age")
head(age.averages[["RNA"]])
head(age.averages[["integrated"]])
df1 <- age.averages[["RNA"]]#
mat <- as.matrix(df1)
dt <- new("ExpressionSet",exprs = mat)
dim(dt)
dt.r <- filter.NA(dt, thres=0.25)
dt.f <- fill.NA(dt.r,mode="mean")
tmp <- filter.std(dt.f,min.std=0.2)
dt.s <- standardise(tmp)
df.s <- dt.s@assayData$exprs
#
m1 <- mestimate(dt.s)
#
set.seed(007)
cl <- mfuzz(dt.s,c=9,m=m1)
#
mfuzz.plot(dt.s,cl,mfrow=c(3,3), new.window= FALSE, time.labels=colnames(dt.s))
#
library(RColorBrewer)
mycol <- c("cyan","yellow","orangered")
mycol <- c("#737fb4","#168c90","#f7c99b","#421c77")
mycolor <- colorRampPalette(mycol)(20)
#
pdf("Others expression changes trends_2(L2-3 IT).pdf", width = 6.5,height = 6)
mfuzz.plot(dt.s,cl,mfrow=c(3,3), new.window= FALSE, time.labels=colnames(dt.s), colo = mycolor)
dev.off()
#
gene_cluster <- data.frame(df.s,cluster=cl$cluster)
write.csv(gene_cluster, 'gene_cluster_Others(L2-3 IT).csv')

#-----------Mfuzz(ET)----------#
library(Mfuzz)
library(Seurat)
library(progeny)
library(tidyr)
library(tibble)
library(dplyr)
library(viridisLite)
library(ggplot2)
#
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
#
scRNA=subset(scRNA,features= rownames(scRNA@assays$RNA@scale.data))
scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
scRNA=ScaleData(scRNA)
#
age.averages <- AverageExpression(scRNA,group.by = "age")
head(age.averages[["RNA"]])
head(age.averages[["integrated"]])
df1 <- age.averages[["RNA"]]#
mat <- as.matrix(df1)
dt <- new("ExpressionSet",exprs = mat)
dim(dt)
dt.r <- filter.NA(dt, thres=0.25)
dt.f <- fill.NA(dt.r,mode="mean")
tmp <- filter.std(dt.f,min.std=0.2)
dt.s <- standardise(tmp)
df.s <- dt.s@assayData$exprs
#
m1 <- mestimate(dt.s)
#
set.seed(007)
cl <- mfuzz(dt.s,c=9,m=m1)
#
mfuzz.plot(dt.s,cl,mfrow=c(3,3), new.window= FALSE, time.labels=colnames(dt.s))
#
library(RColorBrewer)
mycol <- c("cyan","yellow","orangered")
mycol <- c("#737fb4","#168c90","#f7c99b","#421c77")
mycolor <- colorRampPalette(mycol)(20)
#
pdf("Others expression changes trends_2(L2-3 IT).pdf", width = 6.5,height = 6)
mfuzz.plot(dt.s,cl,mfrow=c(3,3), new.window= FALSE, time.labels=colnames(dt.s), colo = mycolor)
dev.off()
#
gene_cluster <- data.frame(df.s,cluster=cl$cluster)
write.csv(gene_cluster, 'gene_cluster_Others(L2-3 ET).csv')

#-----------Mfuzz(InN)----------#
library(Mfuzz)
library(Seurat)
library(progeny)
library(tidyr)
library(tibble)
library(dplyr)
library(viridisLite)
library(ggplot2)
#
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
#
scRNA=subset(scRNA,features= rownames(scRNA@assays$RNA@scale.data))
scRNA@meta.data[["age"]]<-factor(scRNA@meta.data[["age"]], levels=c("Sham","Sham3d","Sham7d","1h","4h","12h","3d","7d"))
#
scRNA=ScaleData(scRNA)
#
age.averages <- AverageExpression(scRNA,group.by = "age")
head(age.averages[["RNA"]])
head(age.averages[["integrated"]])
df1 <- age.averages[["RNA"]]#
mat <- as.matrix(df1)
dt <- new("ExpressionSet",exprs = mat)
dim(dt)
dt.r <- filter.NA(dt, thres=0.25)
dt.f <- fill.NA(dt.r,mode="mean")
tmp <- filter.std(dt.f,min.std=0.2)
dt.s <- standardise(tmp)
df.s <- dt.s@assayData$exprs
#
m1 <- mestimate(dt.s)
#
set.seed(007)
cl <- mfuzz(dt.s,c=9,m=m1)
#
mfuzz.plot(dt.s,cl,mfrow=c(3,3), new.window= FALSE, time.labels=colnames(dt.s))
#
library(RColorBrewer)
mycol <- c("cyan","yellow","orangered")
mycol <- c("#737fb4","#168c90","#f7c99b","#421c77")
mycolor <- colorRampPalette(mycol)(20)
#
pdf("Others expression changes trends_2(L2-3 IT).pdf", width = 6.5,height = 6)
mfuzz.plot(dt.s,cl,mfrow=c(3,3), new.window= FALSE, time.labels=colnames(dt.s), colo = mycolor)
dev.off()
#
gene_cluster <- data.frame(df.s,cluster=cl$cluster)
write.csv(gene_cluster, 'gene_cluster_Others(InN).csv')

#-----------Mfuzz genes set score(IT)----------#
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix)
library(readxl)
library(ggpubr)
#
HIE_gene<-read_xlsx("Mfuzz genes set score.xlsx",col_names = 'gene')
gene<-as.list(HIE_gene)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "HIE")
colnames(AB@meta.data)[11]<-"HIE_Score"
#####IT
AB@meta.data[["celltype"]]<-factor(AB@meta.data[["celltype"]], levels=c("L2-3 IT",  "L4 IT",  "L5 IT",  "L6 IT",  "L6 CLA"))

my_comparisons <- list(c('0','5'),c('1','2'))
my9color <- c('#e4c2ce',  '#EF99C6', '#C25A7F', '#C68778', '#CBAA4D')
library(ggpubr)
p.AddModuleScore <- ggviolin(AB@meta.data, x = "celltype", y = "HIE_Score",
                             color = "celltype",add = 'mean_sd',fill = 'celltype',
                             add.params = list(color = "black")) + 
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 45,hjust = 1)) + 
  NoLegend() + labs(x = '')
p.AddModuleScore
ggsave('Mfuzz genes(IT)_vlo.pdf',width = 4.6,height = 4)

#-----------Mfuzz genes set score(ET)----------#
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix)
library(readxl)
library(ggpubr)
#
HIE_gene<-read_xlsx("Mfuzz genes set score.xlsx",col_names = 'gene')
gene<-as.list(HIE_gene)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "HIE")
colnames(AB@meta.data)[11]<-"HIE_Score"
#####IT
AB@meta.data[["celltype"]]<-factor(AB@meta.data[["celltype"]], levels=c("L2-3 IT",  "L4 IT",  "L5 IT",  "L6 IT",  "L6 CLA"))

my_comparisons <- list(c('0','5'),c('1','2'))
my9color <- c('#FC8A6B',  '#C68778', '#C2A3DA', '#7B6FAB')
library(ggpubr)
p.AddModuleScore <- ggviolin(AB@meta.data, x = "celltype", y = "HIE_Score",
                             color = "celltype",add = 'mean_sd',fill = 'celltype',
                             add.params = list(color = "black")) + 
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 45,hjust = 1)) + 
  NoLegend() + labs(x = '')
p.AddModuleScore
ggsave('Mfuzz genes(ET)_vlo.pdf',width = 4.6,height = 4)

#-----------Mfuzz genes set score(InN)----------#
library(Seurat)
library(tidyverse)
library(cowplot)
library(Matrix)
library(readxl)
library(ggpubr)
#
HIE_gene<-read_xlsx("Mfuzz genes set score.xlsx",col_names = 'gene')
gene<-as.list(HIE_gene)
AB<-AddModuleScore(scRNA, features = gene, ctrl = 100, name = "HIE")
colnames(AB@meta.data)[11]<-"HIE_Score"
#####IT
AB@meta.data[["celltype"]]<-factor(AB@meta.data[["celltype"]], levels=c("L2-3 IT",  "L4 IT",  "L5 IT",  "L6 IT",  "L6 CLA"))

my_comparisons <- list(c('0','5'),c('1','2'))
my9color <- c('#FC8A6B',  '#C68778', '#C2A3DA', '#7B6FAB')
library(ggpubr)
p.AddModuleScore <- ggviolin(AB@meta.data, x = "celltype", y = "HIE_Score",
                             color = "celltype",add = 'mean_sd',fill = 'celltype',
                             add.params = list(color = "black")) + 
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
  scale_color_manual(values = my9color) + 
  scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 45,hjust = 1)) + 
  NoLegend() + labs(x = '')
p.AddModuleScore
ggsave('Mfuzz genes(InN)_vlo.pdf',width = 4.6,height = 4)

#-----------Correlation----------#
###Correlation(ExN)
ExN <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("L6b CT",  "L6 CT",  "L6 CLA", "L5 ET",  "L5 NP",  "L6 IT",  "L5 IT",  "L4 IT",  "L2-3 IT")]
scRNA <- ExN
table(scRNA$age)  
av<-AverageExpression(scRNA,group.by = "celltype", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = "spearman"))
write.csv(cor(av[cg,],method = "spearman"),"cor_ExNcelltype.csv") 
pdf("cor_celltype_ExN.pdf", width = 10,height = 9.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 15,treeheight_col = 15,fontsize = 15)
dev.off()

###Correlation(InN)
InN <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("CR", "Meis2", "Vip",  "Sv2c", "Sst",  "Mme")]
scRNA <- InN
table(scRNA$age)  
av<-AverageExpression(scRNA,group.by = "celltype", assays = "RNA")
av=av[[1]]
head(av)
cg=names(tail(sort(apply(av,1,sd)),1000))
View(av[cg,])
View(cor(av[cg,],method = "spearman"))
write.csv(cor(av[cg,],method = "spearman"),"cor_InNcelltype.csv") 
pdf("cor_celltype_InN.pdf", width = 10,height = 9.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 15,treeheight_col = 15,fontsize = 15)
dev.off()



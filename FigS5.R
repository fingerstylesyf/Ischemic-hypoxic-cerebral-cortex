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

#---------heatmap
scRNA <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("L6b CT",  "L6 CT",  "L5 NP",  "L5 ET",  "L6 CLA", "L6 IT",  "L5 IT",  "L4 IT", "L2-3 IT")]
jjplot <- jjDotPlot(object = scRNA,
                    gene = c("Rfx3","Cux1","Dock5","Lhx8","Rora","Slc30a3","Esrrg","Foxo1","Ntng1","Trps1","Osr1","Nr4a3","Adamts3","Ca3","Rspo1","Synpr","Sulf1","Deptor","Stac","Etv1","Htr2c","Islr2","Tbr1","Syt6","Cplx3","Cobll1"),
                    id = 'celltype',
                    #split.by = 'age',
                    split.by.aesGroup = T,dot.col = c("#2763A2",'white',"#B23030"),
                    point.geom = F,
                    tile.geom = T,ytree = F,cluster.order = c("L6b CT", "L6 CT",  "L5 NP",  "L5 ET",  "L6 CLA", "L6 IT",  "L5 IT",  "L4 IT", "L2-3 IT"),
                    col.min = -1,col.max = 2)#+ coord_flip()
jjplot
ggsave(filename = "heatmap-ExN.pdf", plot = jjplot, device = 'pdf', width = 30, height = 15, units = 'cm')
###
scRNA <- scRNA[,scRNA@meta.data[["celltype"]] %in% c("CR",  "Meis2",  "Sv2c", "Vip",  "Sst",  "Mme")]
jjplot <- jjDotPlot(object = scRNA,
                    gene = c("Tox3","Ptprm","Maf","Npy","Eya4","Mafb","Mkx","Calb1","Nos1","Cdh13","Grin3a","Elfn1","Calb2","Cck","Egfr","Cnr1","Htr3a","Igsf11","Alk","Lamp5","Pdlim5","Chd7","Dcx",'Pbx1',"Prokr2","Meis1","Dcn","Nnat","Rtn1","Nrg1","Tnnt1","Dach1","Plxnd1","Slco2a1","Cxcr4","Ndnf","Nxph4","Kcnk2"),
                    id = 'celltype',
                    #split.by = 'age',
                    split.by.aesGroup = T,dot.col = c("#2763A2",'white',"#B23030"),
                    point.geom = F,
                    tile.geom = T,ytree = F,cluster.order = c("CR", "Meis2",  "Sv2c", "Vip",  "Sst",  "Mme"),
                    col.min = -1,col.max = 2)#+ coord_flip()
jjplot
ggsave(filename = "heatmap-InN.pdf", plot = jjplot, device = 'pdf', width = 40, height = 20, units = 'cm')


#----------highlight
p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='CR'),
                L23_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents="L2-3 IT"),
                L4_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='L4 IT'),
                L5_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='L5 IT'),
                L6_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='L6 IT'),
                L6_CLA=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='L6 CLA'),
                L5_ET=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='L5 ET'),
                L5_NP=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='L5 NP'),
                L6_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='L6 CT'),
                L6b_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='L6b CT'),
                Meis2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='Meis2'),
                Mme=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='Mme'),
                Sst=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='Sst'),
                Sv2c=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='Sv2c'),
                Vip=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='Vip')
                
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#6c584c','#6d9150', '#93B727', '#b4b47b', "#FFE27E",'#7B6FAB','#ff9a94','#C2A3DA',"#CBAA4D","#C68778",'#C25A7F','#FC8A6B','#EF99C6', '#FFE4F5','#6E90D0') 
              
)
p1
ggsave(filename = "Sham.png", plot = p1, device = 'png', width = 18, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='CR'),
                L23_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents="L2-3 IT"),
                L4_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='L4 IT'),
                L5_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='L5 IT'),
                L6_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='L6 IT'),
                L6_CLA=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='L6 CLA'),
                L5_ET=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='L5 ET'),
                L5_NP=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='L5 NP'),
                L6_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='L6 CT'),
                L6b_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='L6b CT'),
                Meis2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='Meis2'),
                Mme=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='Mme'),
                Sst=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='Sst'),
                Sv2c=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='Sv2c'),
                Vip=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='Vip')
                
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#6c584c','#6d9150', '#93B727', '#b4b47b', "#FFE27E",'#7B6FAB','#ff9a94','#C2A3DA',"#CBAA4D","#C68778",'#C25A7F','#FC8A6B','#EF99C6', '#FFE4F5','#6E90D0') 
              
)
p1
ggsave(filename = "Sham3d.png", plot = p1, device = 'png', width = 18, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='CR'),
                L23_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents="L2-3 IT"),
                L4_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='L4 IT'),
                L5_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='L5 IT'),
                L6_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='L6 IT'),
                L6_CLA=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='L6 CLA'),
                L5_ET=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='L5 ET'),
                L5_NP=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='L5 NP'),
                L6_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='L6 CT'),
                L6b_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='L6b CT'),
                Meis2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='Meis2'),
                Mme=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='Mme'),
                Sst=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='Sst'),
                Sv2c=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='Sv2c'),
                Vip=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='Vip')
                
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#6c584c','#6d9150', '#93B727', '#b4b47b', "#FFE27E",'#7B6FAB','#ff9a94','#C2A3DA',"#CBAA4D","#C68778",'#C25A7F','#FC8A6B','#EF99C6', '#FFE4F5','#6E90D0') 
              
)
p1
ggsave(filename = "Sham7d.png", plot = p1, device = 'png', width = 18, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='CR'),
                L23_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents="L2-3 IT"),
                L4_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='L4 IT'),
                L5_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='L5 IT'),
                L6_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='L6 IT'),
                L6_CLA=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='L6 CLA'),
                L5_ET=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='L5 ET'),
                L5_NP=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='L5 NP'),
                L6_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='L6 CT'),
                L6b_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='L6b CT'),
                Meis2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='Meis2'),
                Mme=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='Mme'),
                Sst=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='Sst'),
                Sv2c=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='Sv2c'),
                Vip=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='Vip')
                
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#6c584c','#6d9150', '#93B727', '#b4b47b', "#FFE27E",'#7B6FAB','#ff9a94','#C2A3DA',"#CBAA4D","#C68778",'#C25A7F','#FC8A6B','#EF99C6', '#FFE4F5','#6E90D0') 
              
)
p1
ggsave(filename = "1h.png", plot = p1, device = 'png', width = 18, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='CR'),
                L23_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents="L2-3 IT"),
                L4_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='L4 IT'),
                L5_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='L5 IT'),
                L6_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='L6 IT'),
                L6_CLA=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='L6 CLA'),
                L5_ET=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='L5 ET'),
                L5_NP=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='L5 NP'),
                L6_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='L6 CT'),
                L6b_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='L6b CT'),
                Meis2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='Meis2'),
                Mme=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='Mme'),
                Sst=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='Sst'),
                Sv2c=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='Sv2c'),
                Vip=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='Vip')
                
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#6c584c','#6d9150', '#93B727', '#b4b47b', "#FFE27E",'#7B6FAB','#ff9a94','#C2A3DA',"#CBAA4D","#C68778",'#C25A7F','#FC8A6B','#EF99C6', '#FFE4F5','#6E90D0') 
              
)
p1
ggsave(filename = "4h.png", plot = p1, device = 'png', width = 18, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='CR'),
                L23_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents="L2-3 IT"),
                L4_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='L4 IT'),
                L5_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='L5 IT'),
                L6_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='L6 IT'),
                L6_CLA=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='L6 CLA'),
                L5_ET=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='L5 ET'),
                L5_NP=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='L5 NP'),
                L6_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='L6 CT'),
                L6b_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='L6b CT'),
                Meis2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='Meis2'),
                Mme=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='Mme'),
                Sst=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='Sst'),
                Sv2c=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='Sv2c'),
                Vip=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='Vip')
                
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#6c584c','#6d9150', '#93B727', '#b4b47b', "#FFE27E",'#7B6FAB','#ff9a94','#C2A3DA',"#CBAA4D","#C68778",'#C25A7F','#FC8A6B','#EF99C6', '#FFE4F5','#6E90D0') 
              
)
p1
ggsave(filename = "12h.png", plot = p1, device = 'png', width = 18, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='CR'),
                L23_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents="L2-3 IT"),
                L4_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='L4 IT'),
                L5_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='L5 IT'),
                L6_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='L6 IT'),
                L6_CLA=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='L6 CLA'),
                L5_ET=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='L5 ET'),
                L5_NP=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='L5 NP'),
                L6_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='L6 CT'),
                L6b_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='L6b CT'),
                Meis2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='Meis2'),
                Mme=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='Mme'),
                Sst=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='Sst'),
                Sv2c=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='Sv2c'),
                Vip=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='Vip')
                
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#6c584c','#6d9150', '#93B727', '#b4b47b', "#FFE27E",'#7B6FAB','#ff9a94','#C2A3DA',"#CBAA4D","#C68778",'#C25A7F','#FC8A6B','#EF99C6', '#FFE4F5','#6E90D0') 
              
)
p1
ggsave(filename = "3d.png", plot = p1, device = 'png', width = 18, height = 14.5, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=0.6,cols = "grey85",
              cells.highlight=list(
                CR=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='CR'),
                L23_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents="L2-3 IT"),
                L4_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='L4 IT'),
                L5_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='L5 IT'),
                L6_IT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='L6 IT'),
                L6_CLA=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='L6 CLA'),
                L5_ET=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='L5 ET'),
                L5_NP=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='L5 NP'),
                L6_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='L6 CT'),
                L6b_CT=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='L6b CT'),
                Meis2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='Meis2'),
                Mme=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='Mme'),
                Sst=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='Sst'),
                Sv2c=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='Sv2c'),
                Vip=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='Vip')
                
              ),
              sizes.highlight = c(0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c('#6c584c','#6d9150', '#93B727', '#b4b47b', "#FFE27E",'#7B6FAB','#ff9a94','#C2A3DA',"#CBAA4D","#C68778",'#C25A7F','#FC8A6B','#EF99C6', '#FFE4F5','#6E90D0') 
              
)
p1
ggsave(filename = "7d.png", plot = p1, device = 'png', width = 18, height = 14.5, units = 'cm')


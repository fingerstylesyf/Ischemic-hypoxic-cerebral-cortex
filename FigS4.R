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

####highlight
p1 <- DimPlot(AB, reduction = "umap", pt.size=1.0,cols = "grey85",
              cells.highlight=list(
                MG1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='MG1'),
                #MG2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents="MG2"),
                MG3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='MG3'),
                #MG4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='MG4'),
                #MG5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='MG5'),
                MG6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='MG6'),
                MG7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='MG7'),
                MG8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham")],idents='MG8')
                
              ),
              sizes.highlight = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#2F4858","#ccaf74",'#63b9bd',"#8572ad","#f26f66")
              
)
p1
ggsave(filename = "Sham.png", plot = p1, device = 'png', width = 13, height = 9.8, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=1.0,cols = "grey85",
              cells.highlight=list(
                MG1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='MG1'),
                MG2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents="MG2"),
                MG3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='MG3'),
                MG4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='MG4'),
                MG5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='MG5'),
                MG6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='MG6'),
                MG7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='MG7'),
                MG8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham3d")],idents='MG8')
                
              ),
              sizes.highlight = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#2F4858","#ccaf74",'#63b9bd',"#539075",'#3082b4',"#8572ad","#4962a1","#f26f66")
              
)
p1
ggsave(filename = "Sham3d.png", plot = p1, device = 'png', width = 13, height = 9.8, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=1.0,cols = "grey85",
              cells.highlight=list(
                MG1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='MG1'),
                MG2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents="MG2"),
                MG3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='MG3'),
                MG4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='MG4'),
                MG5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='MG5'),
                MG6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='MG6'),
                MG7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='MG7'),
                MG8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("Sham7d")],idents='MG8')
                
              ),
              sizes.highlight = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#2F4858","#ccaf74",'#63b9bd',"#539075",'#3082b4',"#8572ad","#4962a1","#f26f66")
              
)
p1
ggsave(filename = "Sham7d.png", plot = p1, device = 'png', width = 13, height = 9.8, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=1.0,cols = "grey85",
              cells.highlight=list(
                MG1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='MG1'),
                MG2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents="MG2"),
                #MG3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='MG3'),
                #MG4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='MG4'),
                #MG5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='MG5'),
                #MG6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='MG6'),
                MG7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='MG7'),
                MG8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("1h")],idents='MG8')
                
              ),
              sizes.highlight = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#2F4858","#ccaf74","#4962a1","#f26f66")
              
)
p1
ggsave(filename = "1h.png", plot = p1, device = 'png', width = 13, height = 9.8, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=1.0,cols = "grey85",
              cells.highlight=list(
                MG1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='MG1'),
                MG2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents="MG2"),
                #MG3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='MG3'),
                MG4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='MG4'),
                MG5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='MG5'),
                MG6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='MG6'),
                MG7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='MG7'),
                MG8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("4h")],idents='MG8')
                
              ),
              sizes.highlight = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#2F4858","#ccaf74",'#63b9bd',"#539075",'#3082b4',"#4962a1","#f26f66")
              
)
p1
ggsave(filename = "4h.png", plot = p1, device = 'png', width = 13, height = 9.8, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=1.0,cols = "grey85",
              cells.highlight=list(
                MG1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='MG1'),
                MG2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents="MG2"),
                MG3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='MG3'),
                MG4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='MG4'),
                #MG5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='MG5'),
                #MG6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='MG6'),
                MG7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='MG7'),
                MG8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("12h")],idents='MG8')
                
              ),
              sizes.highlight = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#2F4858","#ccaf74",'#3082b4',"#8572ad","#4962a1","#f26f66")
              
)
p1
ggsave(filename = "12h.png", plot = p1, device = 'png', width = 13, height = 9.8, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=1.0,cols = "grey85",
              cells.highlight=list(
                MG1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='MG1'),
                MG2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents="MG2"),
                MG3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='MG3'),
                MG4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='MG4'),
                MG5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='MG5'),
                MG6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='MG6'),
                MG7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='MG7'),
                MG8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("3d")],idents='MG8')
                
              ),
              sizes.highlight = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#2F4858","#ccaf74",'#63b9bd',"#539075",'#3082b4',"#8572ad","#4962a1","#f26f66")
              
)
p1
ggsave(filename = "3d.png", plot = p1, device = 'png', width = 13, height = 9.8, units = 'cm')

p1 <- DimPlot(AB, reduction = "umap", pt.size=1.0,cols = "grey85",
              cells.highlight=list(
                MG1=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='MG1'),
                MG2=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents="MG2"),
                MG3=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='MG3'),
                MG4=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='MG4'),
                MG5=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='MG5'),
                MG6=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='MG6'),
                MG7=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='MG7'),
                MG8=WhichCells(AB[,AB@meta.data[["age"]] %in% c("7d")],idents='MG8')
                
              ),
              sizes.highlight = c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0), 
              #cols.highlight=list(c5="red", c8="orange")
              cols.highlight=c("#2F4858","#ccaf74",'#63b9bd',"#539075",'#3082b4',"#8572ad","#4962a1","#f26f66")
              
)
p1
ggsave(filename = "7d.png", plot = p1, device = 'png', width = 13, height = 9.8, units = 'cm')

#heatmap
jjplot <- jjDotPlot(object = scRNA,
                    gene = c("Sparc","Nfia","Disc1","Arhgap45","Srgap2","Ets1","Csf1r","Siglec5","Ifngr1","Inpp5d","Entpd1","Havcr2","C3","Tyrobp","Cst3","Snca","Bcl6","Nkain2","Mif","Dpp10","B2m","Fth1","Jund","Aif1","Colec12","F13a1","Fos","Ms4a7","Ppp3ca","Siglec1","Max","Mafb","Cfh","Tnf","Lgmn","Ptprc","Itgax","Ccl3","Ccl4","Hif1a","Il1a","Cd44","Nfkb1","Lgals3","Stat3","Lpl","Itgb2","Olr1","Cd9","Cd63","Plin2","Myo1e","Cdk6","Lpp","Icam1","Nfil3","Plau"),
                    id = 'celltype',
                    #split.by = 'age',
                    split.by.aesGroup = T,dot.col = c("#2763A2",'white',"#B23030"),
                    point.geom = F,
                    tile.geom = T,ytree = F,cluster.order = c("MG1","MG2","MG3","MG4","MG5","MG6","MG7","MG8"),
                    col.min = -1,col.max = 2)#+ coord_flip()

jjplot
ggsave(filename = "heatmap-MG.pdf", plot = jjplot, device = 'pdf', width = 35, height = 10, units = 'cm')












